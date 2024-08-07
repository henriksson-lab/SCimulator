/// 
/// Simulator of single-cell data for scMetaG
/// SCimulator
/// 

use itertools::Itertools;
use rand::distributions::Distribution;
use rand::Rng;
use rand_distr::{WeightedAliasIndex, Uniform, Poisson}; //LogNormal, Normal, 

use ropey::Rope;

type SequenceProvider = dyn Fn() -> Option<Rope>; 
type SequenceProviderFactory = Box<dyn Fn() -> Box<SequenceProvider>>; 

use bio::io::fasta;
use std::path::Path;
use std::io::Error as IOError;
use csv::ReaderBuilder;
use std::error::Error as StdError;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




/// Pick a random genome i.e. simulate a cell;
pub fn new_cell_factory(genomes:Vec<Rope>, prob:Vec<f64>) -> SequenceProviderFactory {
    let dist = WeightedAliasIndex::new(prob).unwrap();
    return Box::new(move || {

        //Pick one genome randomly
        let mut rng = rand::thread_rng();
        let i = dist.sample(&mut rng);
        let this_genome = genomes[i].clone();

        return Box::new(move || {

            //TODO random subsetting
            let r = this_genome.clone();
            return Some(r);
        });
    });
}



/// Pick a random genome i.e. simulate a cell;
/// Support for multiple chromosomes, probability according to size
pub fn new_multichromosome_factory(genomes:Vec<Vec<Rope>>, prob:Vec<f64>) -> SequenceProviderFactory {
    let dist = WeightedAliasIndex::new(prob).unwrap();
    return Box::new(move || {

        //Pick one genome randomly
        let mut rng = rand::thread_rng();
        let i = dist.sample(&mut rng);
        let this_genome = &genomes[i];

        //Pick one chromosome according to size
        let chrom_prob = this_genome.iter().map(|x| x.len_bytes() as f64).collect_vec();
        let chrom_dist = WeightedAliasIndex::new(chrom_prob).unwrap();
        let j = chrom_dist.sample(&mut rng);
        let this_chrom = this_genome[j].clone();

        return Box::new(move || {

            //Get a random subset
            //TODO support wrapping, circular DNA
            let mut rng = rand::thread_rng();
            let from_pos = rng.gen_range(0..this_chrom.len_bytes());

            let sub = Rope::from(this_chrom.byte_slice(from_pos..this_chrom.len_bytes()));
            return Some(sub);
        });
    });
}



/// Pick a random genome proportional to the sizes of the genomes; used to combine multiple
/// chromosomes or plasmids together, assuming equal copy number
pub fn pick_genome_relative_to_size_factory(genomes:Vec<Rope>) -> SequenceProviderFactory {
    let probs = genomes.iter().map(|x| x.len_bytes() as f64).collect_vec();
    return new_cell_factory(
        genomes, 
        probs
    );
}


/// Pick a random genome, equal probability; as in selecting a random cell genome
pub fn pick_genome_equal_probability_factory(genomes:Vec<Rope>) -> SequenceProviderFactory {
    let probs = vec![1.0;genomes.len()];
    return new_cell_factory(
        genomes, 
        probs
    );
}




/// Poisson loading of droplets, having 0+ cells according to lambda
pub fn new_poisson_cell_factory(cell_factory:SequenceProviderFactory, lambda:f64) -> SequenceProviderFactory {
    let dist = Poisson::new(lambda).unwrap();
    return Box::new(move || {

        //Pick N cells, where N is poisson distributed
        let mut rng = rand::thread_rng();
        let num_cell = dist.sample(&mut rng) as u64; //do we really need to round?
        let mut cells = Vec::new();
        for _i in 0..num_cell {
            cells.push(cell_factory());
        }

        return Box::new(move || {
            if cells.len()==0 {
                //If there are no cells, give up sampling
                println!("no cells");
                return None;
            } else {
                //Pick from a random cell. TODO should be relative to cell genome content of each cell
                let mut rng = rand::thread_rng();
                let for_cell = rng.gen_range(0..cells.len());
                let r = cells[for_cell]().clone();
                return r;
            }
       });
   });
}


/// Barcoding and fragmentation of libraries, as performed for Atrandi SPCs
pub fn new_atrandi_barcoding_factory(spc_factory:SequenceProviderFactory, atrandi_barcodes:AtrandiBarcodes) -> SequenceProviderFactory {

    return Box::new(move || {
        //Pick a random barcode
        let combination_bc = atrandi_barcodes.generate_random_bc();

        //Sometimes the first rounds of barcoding do not happen; ~1% of cases. Ignored here. TODO should save ground truth as well
        let bc = Rope::from(format!("{}AGGA{}ACTC{}AAGG{}T", 
            combination_bc[0],
            combination_bc[1],
            combination_bc[2],
            combination_bc[3]
        ));


        let one_spc = spc_factory();
        return Box::new(move || {

            let previous_fragment = one_spc()?;

            //Add barcode. Should in principle add the same barcode to both sides, but likely uncommon to see the right side. We skip this to save performance
            let mut fragment_with_bc=bc.clone();
            fragment_with_bc.append(previous_fragment);
                        
            //NEB fragmentase cuts the fragment in two. Can even happen inside the barcodes!
            let mut rng = rand::thread_rng();        
            let digest_len = rng.gen_range(0..fragment_with_bc.len_bytes()); //replace with better approx
            let mut fragment_digested = Rope::from(fragment_with_bc.byte_slice(0..digest_len));

            //New adapter is added opposite of remaining adapter.
            //the adapter is /5Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTG*T
            fragment_digested.append(Rope::from("GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"));
            return Some(fragment_digested);
       });
   });
}



/// This function returns an atrandi library i.e. a library of cells
pub fn new_atrandi_library_factory(spc_factory:SequenceProviderFactory, num_cells:usize) -> SequenceProviderFactory {
    return Box::new(move || {
        //Generate SPCs
        let mut list_spc = Vec::new();
        for _i in 0..num_cells {
            list_spc.push(spc_factory());
        }

        //For now, same depth TODO
        let prob = list_spc.iter().map(|_x| 1.0).collect_vec();
        let dist = WeightedAliasIndex::new(prob).unwrap();

        return Box::new(move || {
            //Pick one SPC randomly
            let mut rng = rand::thread_rng();
            let i = dist.sample(&mut rng);
            return list_spc[i]();
       });
   });
}




/// This function returns an atrandi library i.e. a library of cells
pub fn new_size_selection_factory(library_factory:SequenceProviderFactory, lower:Option<usize>, upper:Option<usize>) -> SequenceProviderFactory {
    return Box::new(move || {
        let lib = library_factory();

        return Box::new(move || {

            let frag = lib()?;
            let frag_len = frag.len_bytes();

            //Test lower cutoff
            let frag = match lower {
                Some(lower) => {
                    if frag_len >= lower {
                         Some(frag)
                    } else {
                        println!("Failed lower cutoff"); None
                    }
                },
                None => { Some(frag) }
            }?;

            //Test upper cutoff
            let frag = match upper {
                Some(upper) => {
                    if frag_len <= upper {
                            Some(frag)
                    } else {
                        println!("Failed upper cutoff"); None
                    }
                },
                None => { Some(frag) }
            };
            
            return frag;
       });
   });
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Sequencer ///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/// Generate paired libraries like an Illumina sequencer
/// TODO: read errors!
pub fn illumina_generate_read(library:&Box<SequenceProvider>, len_r1:usize, len_r2:usize) -> (String,String) {
    let one_fragment = library();
    match one_fragment {
        Some(one_fragment) => {

            let frag_len=one_fragment.len_bytes();

            //TODO: if fragments shorter than read length, append bogus

            let r1=one_fragment.byte_slice(0..len_r1);
            let r2=one_fragment.byte_slice((frag_len-len_r2)..frag_len);

            //TODO only keep part
            return (r1.as_str().unwrap().into(), r2.as_str().unwrap().into());        
        },
        None => {
            //Try again if we did not get a read
            return illumina_generate_read(library, len_r1, len_r2);
        }
    }
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Helper functions ////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// Read all genomes from a FASTA file, ignoring name
pub fn read_all_seq_from_fasta(f:&Path) -> Result<Vec<Rope>,IOError> {
    let reader = fasta::Reader::from_file(f).expect("Unable to open"); //handle better?
    let mut gens: Vec<Rope> = vec![];
    for result in reader.records() {
        let result_data = &result.unwrap();

        let seq = Rope::from(String::from_utf8_lossy(result_data.seq()));
        //let seq = Rope::from(str::from_utf8(result_data.seq().unwrap()));
        
        println!("{:?}",seq.len_bytes());

        gens.push(seq);
    }
    return Ok(gens);
}


pub fn generate_dna_sequence(length: u16) -> String {
    // generate rand DNA string of length n
    let mut rng = rand::thread_rng();
    let alphabet: Vec<&str> = vec!["A", "C", "G", "T"];
    let between = Uniform::from(0..alphabet.len());

    let mut seq: String = String::new();
    for _ in 0..length {
        let s = alphabet[between.sample(&mut rng)];
        seq += s;
    }
    seq
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// List of barcodes for one round
pub struct BarcodeWhitelist {
    list: Vec<Rope>    //List for alignment; not sure if worth having separate from set
}

/// Structure for Atrandi combinatorial barcodes
pub struct AtrandiBarcodes {
    rounds: Vec<BarcodeWhitelist>
}


impl AtrandiBarcodes {

    /// Read dictionary of Atrandi barcodes from file
    pub fn read_atrandi_barcodes(filename:&str) -> Result<AtrandiBarcodes, Box<dyn StdError>> {
        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(filename)?;
        let mut bcs_for_well = vec![vec![] as Vec<String>; 4];
        for result in rdr.records() {
            let record = result?;
            let pos=&record[0];
            let bc=&record[2];
            let pos_int = pos.parse::<usize>().unwrap() - 1;
            bcs_for_well[pos_int].push(String::from(bc));        
        }

        let whitelists = bcs_for_well.iter().map(|w| BarcodeWhitelist {
            list: w.iter().map(|x| Rope::from_str(x)).collect_vec()
        }  ).collect();
        
        Ok(AtrandiBarcodes {rounds: whitelists})
    }


    pub fn generate_random_bc(&self) -> Vec<Rope> {

        let mut rng = rand::thread_rng();
        
        let random_bcs = self.rounds.iter().
            map(|round| round.list[rng.gen_range(0..round.list.len())].clone()).
            collect_vec();

        return random_bcs;
    }

}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





fn simulate_atrandi_wgs() {



    //Load all genomes; ignoring for now copy number of plasmids
    let list_fasta = vec![
        Path::new("/Users/mahogny/Desktop/rust/playground/genomes/Bacillus cereus (ATCC 10987) NC_003909.8.fasta")
    ];
    let genomes = list_fasta.iter().map(
        |f| read_all_seq_from_fasta(f).expect("Unable to open")).
        collect_vec();
    
    let genome_probs = vec![1.0; genomes.len()];
    let cell_factory = new_multichromosome_factory(
        genomes, 
        genome_probs
    );


    //An SPC is 0+ cells together
    let spc_factory = new_poisson_cell_factory(
        cell_factory,
        10.0 //0.3
    ); 

    //Adds barcodes and fragments the DNA
    let atrandi_barcodes: AtrandiBarcodes = AtrandiBarcodes::read_atrandi_barcodes("/Users/mahogny/Desktop/rust/playground/atrandi_bc.csv").expect("Unable to open barcode file");
    let barcoding_factory = new_atrandi_barcoding_factory(spc_factory, atrandi_barcodes);
    
    //Makes a library from SPCs of different depth
    let spc_library_factory = new_atrandi_library_factory(
        barcoding_factory, 
        1000  //Number of SPCs; should also take distribution of reads across them, and background; could also be separate
    ); 

    let size_selected_library_factory = new_size_selection_factory(
        spc_library_factory, 
        Some(200), None
    ); 

    

    //Instanciate a library
    let the_library = size_selected_library_factory();

    //Loop to generate enough reads
    for _i in 0..100 {

        let (r1,r2) = illumina_generate_read(
            &the_library,
            150,
            150
        );
    
    
        println!("{}\n{}\n", r1, r2);
        //println!("{}", r1);
    
    }


}



fn main() {

    simulate_atrandi_wgs();

}

