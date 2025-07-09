//use core::num;
use std::{collections::HashMap, fs::File, io::{self, BufRead, Write}, path::Path};
use std:: sync::mpsc::channel;
use clap::{Arg, ArgMatches, Command};
use chrono::{DateTime, Local};
use rayon::prelude::*;
//use rand::Rng;
use gchemol::prelude::*;
use gchemol::Molecule;
//use spdkit::prelude::*;

use rustamath_mnmz::amoeba;
use rest_tensors::{matrix_blas_lapack::_dgemm_full, BasicMatrix, MatrixFull};
use libm::{cos, sin};
use liblbfgs::lbfgs;
use lazy_static::lazy_static;

use regex::Regex;

pub fn parse_input() -> ArgMatches {
    Command::new("similarity_check")
        .version("0.1")                                                                                                                                                      
        .author("Igor Ying Zhang <igor_zhangying@fudan.edu.cn>")
        .about("Similarity check for a batch of molecular geometries")
        .arg(Arg::new("mode")
             .short('m')
             .long("mode")
             .default_value("0")
             .value_name("mode")
             .help("The filter mode:
 0: filter out a set of unique geometries in <input_file> wrt the <similarity_threshold> 
 1: Step 1, locate the geometry with the minal energy <min>. 
    Step 2, compare with <min> and filter out a set of molecules wrt the <similarity_threshold> and <energy_threshold>
 2: same as mode 0 to filter out a set of gemoetries, but exclude the geometries similar to those in the training set <training_file>
 3: same as mode 2, but continude the filter procedue until the amount of selected geometries reaches a given number <number_threhold>
")
            )
        .arg(Arg::new("input_file")
             .short('i')
             .long("input-file")
             .value_name("input_file")
             .default_value("geometries.xyz")
             .help("A batch of molecular geometries stored in XYZ format")
            )
        .arg(Arg::new("output_file")
             .short('o')
             .long("output-file")
             .value_name("output_file")
             .default_value("filter_out.xyz")
             .help("A batch of molecular geometries filtered out with similarity less than threshold")
            )
        .arg(Arg::new("training_file")
             .short('r')
             .long("training_file")
             .value_name("training_file")
             .default_value("training.xyz")
             .help("For mode 2 and 3, it needs a batch of molecular geometries stored in XYZ format, used as a constain")
            )
        .arg(Arg::new("similarity_threshold")
             .short('t')
             .long("similarity_threshold")
             .value_name("threshold")
             .default_value("0.2")
             .help("Threshold for similarity check based on the root mean squared error (RMSE) between geometries")
            )
        .arg(Arg::new("energy_threshold")
             .short('e')
             .long("energy_threshold")
             .default_value("0.1")
             .value_name("energy_threshold")
             .default_value("0.001")
             .help("For mode = 1, filter out the geometries within the energy window: [minimal energy, minimal energy + energy_threshold]")
            )
        .arg(Arg::new("number_output")
             .short('g')
             .long("number_output")
             .default_value("1000")
             .value_name("number_output")
             .help("For mode = 3, filter out a group of geometries with the total number of geometries = <number_output>")
            )
        .arg(Arg::new("block_size")
             .short('s')
             .long("block_size")
             .value_name("block_size")
             .default_value("100")
             .help("Block-based corase grain approximation:\n Filter the molecules block by block and the similarity check is performed within the block only.")
            )
        .arg(Arg::new("displace")
             .short('d')
             .long("displace")
             .default_value("0.001")
             .value_name("displace")
             .help("Displacement to generate numerical gradients used for lbfgs optimizaiton")
            )
        .arg(Arg::new("num_threads")
             .short('n')
             .long("num_threads")
             .default_value("8")
             .value_name("num_threads")
             .help("Setting the number of threads used for rayon parallelism")
            )
        .arg(Arg::new("fine_filter")
             .short('f')
             .long("fine_filter")
             .action(clap::ArgAction::SetTrue)
             .default_value("false")
             .value_name("fine_filter")
             .help("Turn on the fine filter after block-based corase grain filter procedure")
            )
        .get_matches()
}

pub const SPECIES_NAME: [&str; 72] = ["H", "He",
        "Li","Be","B", "C",
        "N", "O", "F", "Ne",
        "Na","Mg","Al","Si",
        "P", "S", "Cl","Ar",
        "K", "Ca","Ga","Ge",
        "As","Se","Br","Kr",
        "Sc","Ti","V", "Cr","Mn",
        "Fe","Co","Ni","Cu","Zn",
        "Rb","Sr","In","Sn",
        "Sb","Te","I", "Xe",
        "Y", "Zr","Nb","Mo","Tc",
        "Ru","Rh","Pd","Ag","Cd",
        "Cs","Ba","Tl","Pb",
        "Bi","Po","At","Rn",
        "La","Hf","Ta","W", "Re",
        "Os","Ir","Pt","Au","Hg"
        ];

pub const MASS_CHARGE: [(f64,f64);72] = [(1.00794,1.0), (4.002602,2.0),
   (6.9410, 3.0),  (9.0122, 4.0), (10.8110, 5.0), (12.0107, 6.0), 
  (14.0067, 7.0), (15.9994, 8.0), (18.9984, 9.0), (20.1797,10.0),
  (22.9897,11.0), (24.3050,12.0), (26.9815,13.0), (28.0855,14.0), 
  (30.9738,15.0), (32.0650,16.0), (35.4530,17.0), (39.9480,18.0),
  (39.0983,19.0), (40.0780,20.0), (69.7230,31.0), (72.6400,32.0), 
  (74.9216,33.0), (78.9600,34.0), (79.9040,35.0), (83.7980,36.0),
  (44.9559,21.0), (47.8670,22.0), (50.9415,23.0), (51.9961,24.0), (54.9380,25.0),
  (55.8450,26.0), (58.9332,27.0), (58.6934,28.0), (63.5460,29.0), (65.3800,30.0),
  (85.4678,37.0), (87.6200,38.0),(114.8180,49.0),(118.7100,50.0),
 (121.7600,51.0),(127.6000,52.0),(126.9045,53.0),(131.2930,54.0),
  (88.9058,39.0), (91.2240,40.0), (92.9064,41.0), (95.9600,42.0), (98.0000,43.0),
 (101.0700,44.0),(102.9055,45.0),(106.4200,46.0),(107.8682,47.0),(112.4110,48.0),
 (132.9054,55.0),(137.3270,56.0),(204.3833,81.0),(207.2000,82.0),
 (208.9804,83.0),(209.0000,84.0),(210.0000,85.0),(222.0000,86.0),
 (138.9055,57.0),(178.4900,72.0),(180.9479,73.0),(183.8400,74.0),(186.2070,75.0),
 (190.2300,76.0),(192.2170,77.0),(195.0840,78.0),(196.9666,79.0),(200.5900,80.0)
   ];

lazy_static!{
    pub static ref SPECIES_INFO: HashMap<&'static str, &'static (f64,f64)> = {
        let mut m = HashMap::new();
        SPECIES_NAME.iter().zip(MASS_CHARGE.iter()).for_each(|(name,info)| {
            m.insert(*name,info);
        });
        m
    };
}

pub fn formated_element_name(elem: &String) -> String {
    // make sure the first letter of the element name is in uppercase.
    let tmp_elem = elem.to_lowercase();
    let mut c = tmp_elem.chars();
    match c.next() {
        None => String::new(),
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
    }
}


pub fn get_mass(elem_list: &Vec<String>) -> Vec<f64> {

    let mut final_list: Vec<f64> = vec![];

    for elem in elem_list {
        let formated_elem = formated_element_name(&elem);
        let tmp_result = SPECIES_INFO.get(formated_elem.as_str());
        if let Some(tmp_value) = tmp_result {
            final_list.push((**tmp_value).0);
        } else {
            panic!("Specify unknown element: {}. Please check the input file", elem);
        }
    }

    final_list
}

pub fn init_timing() -> DateTime<Local> {
    Local::now()
}

pub fn timing(dt0: &DateTime<Local>, iprint: Option<&str>) -> DateTime<Local> {
//pub fn timing(dt0: &DateTime<Local>, iprint: Option<&str>) {
    let dt1 = Local::now();
    match iprint {
        Some(header) => {
            let timecost1 = (dt1.timestamp_millis()-dt0.timestamp_millis()) as f64 /1000.0;
            println!("{:30} cost {:6.2} seconds", header, timecost1);
        },
        None => {()}
    }
    dt1
}

pub fn sort_geometry(xyz: &MatrixFull<f64>, elem: &Vec<String>) -> HashMap<String, MatrixFull<f64>> {
    let mut xyz_hash: HashMap<String, MatrixFull<f64>> = HashMap::new();
    xyz.iter_columns_full().zip(elem.iter()).for_each(|(col, ele)| {
        if xyz_hash.contains_key(ele) {
            let submatrix = xyz_hash.get_mut(ele).unwrap();
            let mut shape: [usize; 2] = submatrix.size().try_into().unwrap();
            shape[1] += 1;
            let mut tmp_vec = submatrix.data.clone();
            for x in col {
                tmp_vec.push(*x);
            };
            *submatrix = MatrixFull::from_vec(shape, tmp_vec).unwrap();
        } else {
            let tmp_vec = col.iter().map(|x| *x).collect::<Vec<f64>>();
            xyz_hash.insert(ele.clone(), MatrixFull::from_vec([3,1], tmp_vec).unwrap());
        }
    });

    xyz_hash
}

pub fn rotation_3d(xyz: &MatrixFull<f64>, rot: &[f64]) -> MatrixFull<f64> {

    let a = rot[0];
    let cos_a = cos(a);
    let sin_a = sin(a);
    let b = rot[1];
    let cos_b = cos(b);
    let sin_b = sin(b);
    let c = rot[2];
    let cos_c = cos(c);
    let sin_c = sin(c);

    let rot_a = MatrixFull::from_vec([3,3], vec![
        1.0, 0.0, 0.0, 0.0, cos_a, -sin_a, 0.0, sin_a, cos_a
    ]).unwrap();
    let rot_b = MatrixFull::from_vec([3,3], vec![
        cos_b, 0.0, sin_b, 0.0, 1.0, 0.0, -sin_b, 0.0, cos_b
    ]).unwrap();
    let rot_c = MatrixFull::from_vec([3,3], vec![
        cos_c, -sin_c, 0.0, sin_c, cos_c, 0.0, 0.0, 0.0, 1.0
    ]).unwrap();

    let mut rot_o = MatrixFull::new([3,3],0.0);
    _dgemm_full(&rot_a, 'N', &rot_b, 'N', &mut rot_o, 1.0, 0.0);
    _dgemm_full(&rot_o.clone(), 'N', &rot_c, 'N', &mut rot_o, 1.0, 0.0);
    //_dgemm_full(&rot_overall, &rot_c, &mut rot_overall, 1.0, 0.0);

    let mut rot_xyz = xyz.clone();
    _dgemm_full(&rot_o, 'N', xyz, 'N', &mut rot_xyz, 1.0, 0.0);

    rot_xyz

}

pub fn diff_mols(mol1: &MatrixFull<f64>, mol2: &MatrixFull<f64>) -> f64 {
    let diff = mol1.iter().zip(mol2.iter()).map(|(a, b)| (a - b).powi(2)).sum::<f64>();

    (diff/(mol1.data().len() as f64)).sqrt()
}

pub fn similarity_check(mol: &MatrixFull<f64>, ref_mol: &MatrixFull<f64>, threshold: f64, atm_type: &[String]) -> (bool, f64) {

    let mut mol1 = matrix2mol(mol, atm_type);
    let mut mol2 = matrix2mol(ref_mol, atm_type);

    mol1.rebond();
    mol2.rebond();


    //if mol1.fingerprint().eq(&mol2.fingerprint()) {
    //    return (true, -1.0)
    //} else {
    //    let loss = mol1.superimpose_onto(&mol2, None);

    //    if loss < threshold {
    //        return (true, loss)
    //    } else {
    //        return (false, loss)
    //    }
    //}
    let loss = mol1.superimpose_onto(&mol2, None);

    if loss < threshold {
        return (true, loss)
    } else {
        return (false, loss)
    }
}

pub fn similarity_optimization_lbfgs(mol: &MatrixFull<f64>, ref_mol: &MatrixFull<f64>, displace: f64, threshold: f64) -> (MatrixFull<f64>, [f64;3], f64, f64) {

    let check = |x: &[f64]| {
        let rot = rotation_3d(mol, x);
        diff_mols(&rot, ref_mol)
    };

    let fmin_org = diff_mols(&mol, &ref_mol);

    if fmin_org <= threshold {
        return (mol.clone(), [0.0;3], fmin_org, fmin_org)
    }

    //let mut rot_batch = [[0.0;3], [1.0;3], [0.5;3], [0.0, 0.5, 1.0]];
    let mut rot_batch = [[0.0;3]];
    let mut fmin_final = fmin_org;
    let mut rot_final = [0.0;3];
    let mut count = 0;

    rot_batch.iter_mut().for_each(|rotation| {
        let return_state = lbfgs().minimize(rotation, 
        |x: &[f64], gx: &mut [f64]| {
                    let rot = rotation_3d(mol, x);
                    let loss = diff_mols(&rot, ref_mol);

                    for i in 0..3 {
                        let mut rot_x = x.to_vec();
                        rot_x[i] += displace;
                        let rot_p = rotation_3d(mol, &rot_x);
                        let loss_p0 = diff_mols(&rot_p, ref_mol);
                        rot_x[i] -= 2.0*displace;
                        let rot_p = rotation_3d(mol, &rot_x);
                        let loss_p1 = diff_mols(&rot_p, ref_mol);
                        gx[i] = 0.5*(loss_p0 - loss_p1) / displace;
                    }

                    Ok(loss)
                },
            |_| {
                //println!("Iteration {}, Evaluation: {}", &prgr.niter, &prgr.neval);
                //println!(" xnorm = {}, gnorm = {}, step = {}",
                //    &prgr.xnorm, &prgr.gnorm, &prgr.step
                //);
                false
            },
        );

        let fmin = if let Err(_) = return_state {
            fmin_org
        } else {
            check(rotation)
        }; 


        if (fmin - fmin_final) < -1.0e-3 {
            if count > 0 {println!("For benchmark: fmin_org: {:16.8}, fmin_prev: {:16.8}, fmin: {:16.8}", fmin_org, fmin_final, fmin)};
            fmin_final = fmin;
            rot_final = rotation.clone();
            count += 1;
        };
    });


    //(rotation_3d(mol, &rot_final), rot_final, fmin_org, fmin_final)
    (mol.clone(), rot_final, fmin_org, fmin_final)
}

pub fn similarity_optimization(mol: &MatrixFull<f64>, ref_mol: &MatrixFull<f64>) -> (MatrixFull<f64>, Vec<f64>) {

    let check = |x: &[f64]| {
        let rot = rotation_3d(mol, x);
        diff_mols(&rot, ref_mol)
    };

    let fmin_org = check(&[0.0, 0.0, 0.0]);

    let initial_guess = [[0.0,0.0,0.0], [1.0, 1.0, 1.0], [0.5, -0.5, 0.5]];
    let mut fmin_o = fmin_org;
    let mut min_o = vec![0.0, 0.0, 0.0];

    initial_guess.iter().for_each(|initial_guess| {
        let (min, fmin, nr_iterations) = amoeba(
            check, initial_guess, 0.05, 0.05, 500
        );
        if fmin < fmin_o {
            fmin_o = fmin;
            min_o = min;
            println!("fmin_org: {:16.8}, fmin: {:16.8}, nr_iterations: {:5}", fmin_org, fmin_o, nr_iterations);
        }
        
    });

    //println!("fmin_org: {:16.8}, fmin: {:16.8}, nr_iterations: {:5}", fmin_org, fmin, nr_iterations);

    (rotation_3d(mol, &min_o),min_o)
}

/// given a group of atoms, the mass of them are stored in mass_list: &Vec<f64>, their positions are in &MatrixFull<f64>, calculate the mass center
pub fn mass_center(atoms: &MatrixFull<f64>, mass_list: &Vec<f64>) -> [f64;3] {

    let mut mass_center = [0.0, 0.0, 0.0];
    let mut total_mass = 0.0;

    atoms.iter_columns_full().zip(mass_list.iter()).for_each(|(atom, mass)| {
        mass_center.iter_mut().zip(atom.iter()).for_each(|(a, b)| *a += *b * mass);
        total_mass += mass;
    });

    mass_center[0] /= total_mass;
    mass_center[1] /= total_mass;
    mass_center[2] /= total_mass;

    mass_center
}

pub fn shift_to_mass_center(mol: &mut MatrixFull<f64>, mass_list: &Vec<f64>) {

    let mass_center = mass_center(mol, mass_list);

    mol.iter_columns_full_mut().for_each(|x| {
        x.iter_mut().zip(mass_center.iter()).for_each(|(a, b)| *a -= *b)
    });
}

//pub fn ordering_check(mol1: &MatrixFull<f64>, atm_type: &[String], mol2: &MatrixFull<f64>, atm_type2: &[String]) -> bool {
//
//    let mut mol1_tmp = mol1.clone();
//    let mut mol2_tmp = mol2.clone();
//
//    shift_to_mass_center(&mut mol1_tmp, &get_mass(&atm_type));
//    shift_to_mass_center(&mut mol2_tmp, &get_mass(&atm_type2));
//
//    false
//} 

pub fn shift_wrt_ref_atom(mol: &mut MatrixFull<f64>, ref_atom: usize) {
    let ref_atom_coords = mol.iter_column(ref_atom).map(|x| *x).collect::<Vec<f64>>();

    mol.iter_columns_full_mut().for_each(|x| {
        x.iter_mut().zip(ref_atom_coords.iter()).for_each(|(a, b)| *a -= *b)
    });
}

pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>> where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

//pub fn import_xyz_in_batch_gchemol(filename: &str) -> Result<Vec<Molecule>{
//    use gchemol::io;
//
//    let mols: Vec<_> = io::read(filename)?.collect()
//
//
//}

pub fn matrix2xyz(matr: &MatrixFull<f64>, atm_type: &[String]) -> String {
    let mut xyz = format!("{}\n\n", atm_type.len());

    matr.iter_columns_full().zip(atm_type.iter()).for_each(|(coord, atm_type)| {
    
        xyz = format!("{} {:6} {:16.8} {:16.8} {:16.8}\n", xyz, atm_type, coord[0], coord[1], coord[2]);

    } );

    xyz
}

pub fn matrix2mol(matr: &MatrixFull<f64>, atm_type: &[String]) -> Molecule {
    let xyz = matrix2xyz(matr, atm_type);

    Molecule::from_str(&xyz, "text/xyz").unwrap()

}

pub fn import_xyz_in_batch(filename: &str) -> anyhow::Result<(Vec<MatrixFull<f64>>, Option<Vec<f64>>, Vec<Vec<String>>, Option<Vec<MatrixFull<f64>>>)> {

    let re_atmnum = Regex::new(r"^\s*(?P<atm_num>\d+)\s*$").unwrap();

    let re_coord = Regex::new(
      r"(?P<elem>\w{1,2})\s+(?P<x>[\+-]?\d+.\d+)\s+(?P<y>[\+-]?\d+.\d+)\s+(?P<z>[\+-]?\d+.\d+)"
    ).unwrap();

    let re_pos_mom_for = Regex::new(
      r"(?P<elem>\w{1,2})\s+(?P<px>[\+-]?\d+.\d+)\s+(?P<py>[\+-]?\d+.\d+)\s+(?P<pz>[\+-]?\d+.\d+)\s+(?P<mx>[\+-]?\d+.\d+)\s+(?P<my>[\+-]?\d+.\d+)\s+(?P<mz>[\+-]?\d+.\d+)\s+(?P<fx>[\+-]?\d+.\d+)\s+(?P<fy>[\+-]?\d+.\d+)\s+(?P<fz>[\+-]?\d+.\d+)"
    ).unwrap();

    let re_energy = Regex::new(
      r"\s+energy=\s*(?P<energy>[\+-]?\d+.\d+)\s*"
    ).unwrap();

    let file = File::open(filename)?;
    let lines = io::BufReader::new(file).lines();

    let mut mol_batch: Vec<MatrixFull<f64>> = vec![];
    let mut mol_energy: Vec<f64> = vec![];
    let mut mol_momfor: Vec<MatrixFull<f64>> = vec![];
    let mut mol_atmtyp: Vec<Vec<String>> = vec![];

    let mut title_line = false;
    let mut atm_num = 0_usize;
    let mut tmp_coord: Vec<f64> = vec![];
    let mut atm_type : Vec<String> = vec![];
    let mut tmp_mf: Vec<f64> = vec![];

    for line in lines {
        if let Ok(ip) = line {
            if title_line {
                //println!("title line: {}", &ip);
                title_line = false;
                tmp_coord = Vec::with_capacity(atm_num*3);
                atm_type = Vec::with_capacity(atm_num);
                tmp_mf = Vec::with_capacity(atm_num*6);
                if re_energy.is_match(&ip) {
                    let caps = re_energy.captures(&ip).unwrap();
                    let energy = caps.name("energy").unwrap().as_str().parse::<f64>().unwrap();
                    //println!("energy: {}", energy);
                    mol_energy.push(energy);
                } else {
                    println!("the energy is not found in the title line");
                }
            } else {
                if re_atmnum.is_match(&ip) {
                    //println!("atm_num line: {}", &ip);
                    if tmp_coord.len()>0 {
                        mol_batch.push(MatrixFull::from_vec([3,atm_num],tmp_coord.clone()).unwrap());
                        mol_atmtyp.push(atm_type.clone());
                        if tmp_mf.len()>0 {
                            mol_momfor.push(MatrixFull::from_vec([6, atm_num], tmp_mf.clone()).unwrap());
                        }
                    };
                    let caps = re_atmnum.captures(&ip).unwrap();
                    atm_num = caps.name("atm_num").unwrap().as_str().parse::<usize>().unwrap();
                    title_line = true;
                } else if re_pos_mom_for.is_match(&ip) {
                    let caps = re_pos_mom_for.captures(&ip).unwrap();
                    //println!("debug: {}", &ip);
                    //println!("debug: {:?}", &caps);
                    let px = caps.name("px").unwrap().as_str().parse::<f64>().unwrap();
                    let py = caps.name("py").unwrap().as_str().parse::<f64>().unwrap();
                    let pz = caps.name("pz").unwrap().as_str().parse::<f64>().unwrap();
                    tmp_coord.extend([px, py, pz]);
                    let elem = caps.name("elem").unwrap().as_str();
                    atm_type.push(elem.to_string());
                    let mx = caps.name("mx").unwrap().as_str().parse::<f64>().unwrap();
                    let my = caps.name("my").unwrap().as_str().parse::<f64>().unwrap();
                    let mz = caps.name("mz").unwrap().as_str().parse::<f64>().unwrap();
                    let fx = caps.name("fx").unwrap().as_str().parse::<f64>().unwrap();
                    let fy = caps.name("fy").unwrap().as_str().parse::<f64>().unwrap();
                    let fz = caps.name("fz").unwrap().as_str().parse::<f64>().unwrap();
                    tmp_mf.extend([mx, my, mz, fx, fy, fz]);
                } else if re_coord.is_match(&ip) {
                    //println!("coord line: {}", &ip);
                    let caps = re_coord.captures(&ip).unwrap();
                    let x = caps.name("x").unwrap().as_str().parse::<f64>().unwrap();
                    let y = caps.name("y").unwrap().as_str().parse::<f64>().unwrap();
                    let z = caps.name("z").unwrap().as_str().parse::<f64>().unwrap();
                    tmp_coord.extend([x, y, z]);
                    let elem = caps.name("elem").unwrap().as_str();
                    atm_type.push(elem.to_string());
                } else {
                    println!("unknown line: {}", &ip);
                }
            }
        }
    }

    if tmp_coord.len()>0 {
        mol_batch.push(MatrixFull::from_vec([3,atm_num],tmp_coord.clone()).unwrap());
        mol_atmtyp.push(atm_type.clone());
        if tmp_mf.len()>0 {
            mol_momfor.push(MatrixFull::from_vec([6, atm_num], tmp_mf.clone()).unwrap());
        }
    };


    //if mol_batch.len() == mol_energy.len() {
    //    if mol_momfor.len() == mol_batch.len() {
    //        Ok((mol_batch, Some(mol_energy), mol_atmtyp, Some(mol_momfor)))
    //    } else {
    //        Ok((mol_batch, Some(mol_energy), mol_atmtyp, None))
    //    }
    //    //Ok((mol_batch, Some(mol_energy), mol_atmtyp))
    //} else {
    //    println!("The number of molecules and energies do not match\n Molecules: {:16.8}; Energies: {:16.8}", mol_batch.len(), mol_energy.len());
    //    Ok((mol_batch, None, mol_atmtyp))
    //}
    let out_energy = if mol_batch.len() == mol_energy.len() {
        Some(mol_energy)
    } else {
        println!("The number of molecules and energies do not match\n Molecules: {:16.8}; Energies: {:16.8}", mol_batch.len(), mol_energy.len());
        None
    };
    let out_monfor = if mol_momfor.len() == mol_batch.len() {
        Some(mol_momfor)
    } else {
        println!("The number of molecules and forces do not match\n Molecules: {:16.8}; Forces: {:16.8}", mol_batch.len(), mol_momfor.len());
        None
    };

    //Ok((mol_batch, mol_energy, atm_type))

    Ok((mol_batch, out_energy, mol_atmtyp, out_monfor))

}

pub fn filter_batch(mols_batch: &[MatrixFull<f64>], _displace: f64, threshold: f64, atm_type: &[String]) -> Vec<MatrixFull<f64>> {
    //let mols_batch_curr = mols_batch.clone();

    let mut selected_mols: Vec<MatrixFull<f64>> = vec![];

    let mut mols_batch_next = mols_batch.to_vec();

    while mols_batch_next.len() > 0 {

        println!("====================================================================");
        println!("    Remaining molecules: {:5}", mols_batch_next.len());
        println!("====================================================================");

        let ref_mol = mols_batch_next.remove(0);
       selected_mols.push(ref_mol.clone()); 

        let mols_batch_curr = mols_batch_next.clone();

        mols_batch_next.clear();

        mols_batch_curr.iter().for_each(|mol| {
            //let (rot_mol,_, fmin_org, fmin) = similarity_optimization_lbfgs(mol, &ref_mol, displace, threshold);
            //if fmin > threshold {
            //    mols_batch_next.push(rot_mol);
            //}
            //println!("fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);
            let (flag, loss) = similarity_check(mol, &ref_mol, threshold, atm_type);
            if ! flag {
                mols_batch_next.push(mol.clone());
                println!("loss: {:16.8}", loss);
            }
        });
    }

    selected_mols

}

pub fn filter_batch_rayon(mols_batch: &[MatrixFull<f64>], _displace: f64, threshold: f64, atm_type: &[String]) -> Vec<MatrixFull<f64>> {

    //let thread_number = rayon::current_num_threads();

    let mut selected_mols: Vec<MatrixFull<f64>> = vec![];

    let mut mols_batch_next = mols_batch.to_vec();

    while mols_batch_next.len() > 1 {

        println!("====================================================================");
        println!(" Selected molecules: {:5}. Remaining molecules: {:5}", selected_mols.len(), mols_batch_next.len());
        println!("====================================================================");

        let ref_mol = mols_batch_next.remove(0);
        selected_mols.push(ref_mol.clone()); 

        let mols_batch_curr = mols_batch_next.clone();
        mols_batch_next.clear();


        let (sender, receiver) = channel();

        let mols_size = mols_batch_curr.len();
        let thrds_num = rayon::current_num_threads();
        let mut  batch_size_thread = mols_size/thrds_num;
        let left_issue = mols_size%thrds_num;
        if left_issue > 0 {
            batch_size_thread += 1;
        }

        mols_batch_curr.par_chunks(batch_size_thread).for_each_with(sender,|s, sub_mols_batch| {
            let mut sub_mols_batch_next = vec![];
            sub_mols_batch.iter().for_each(|mol| {
                //let (rot_mol, _, fmin_org, fmin) = similarity_optimization_lbfgs(mol, &ref_mol, displace, threshold);

                //if fmin > threshold {
                //    sub_mols_batch_next.push(rot_mol);
                //}

                //if let Some(thread_index) = rayon::current_thread_index() {
                //    if thread_index == 0 {
                //        println!("thread 0 => fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);
                //    }
                //} else {
                //    println!("fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);
                //}
                let (flag, loss) = similarity_check(mol, &ref_mol, threshold, atm_type);
                if ! flag {
                    sub_mols_batch_next.push(mol.clone());
                    //println!("fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);
                };
                if let Some(thread_index) = rayon::current_thread_index() {
                    if thread_index == 0 {
                        println!("loss: {:16.8}", loss);
                    }
                } 
            });
            s.send(sub_mols_batch_next).unwrap();
        });

        receiver.iter().for_each(|sub_mols_batch_next| {
            mols_batch_next.extend(sub_mols_batch_next);
        });

        //mols_batch_curr.par_iter().for_each_with(sender, |s, mol| {
        //    let (rot_mot, _, fmin_org, fmin) = similarity_optimization_lbfgs(mol, &ref_mol, displace);
        //    if let Some(thread_index) = rayon::current_thread_index() {
        //        if thread_index == 0 {
        //            println!("fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);
        //        }
        //    };
        //    let out_mol = if fmin > threshold {
        //        Some(rot_mot)
        //    } else {
        //        None
        //    };
        //    s.send(out_mol).unwrap();
        //});
        //receiver.iter().for_each(|out_mol| {
        //    if let Some(mol) = out_mol {
        //        mols_batch_next.push(mol);
        //    }
        //});
    };

    if let Some(mol) = mols_batch_next.pop() {  
        selected_mols.push(mol);
    }

    selected_mols
}

pub fn filter_batch_rayon_with_energy(mols_batch: &[MatrixFull<f64>], energy_batch: &[f64], _displace: f64, threshold: f64, atm_type: &[String]) -> (Vec<MatrixFull<f64>>, Vec<f64>) {

    //let thread_number = rayon::current_num_threads();

    let mut selected_mols: Vec<MatrixFull<f64>> = vec![];
    let mut selected_enys: Vec<f64> = vec![];

    let mut mols_batch_next = mols_batch.to_vec();
    let mut enys_batch_next = energy_batch.to_vec();

    while mols_batch_next.len() > 1 {

        println!("====================================================================");
        println!(" Selected molecules: {:5}. Remaining molecules: {:5}", selected_mols.len(), mols_batch_next.len());
        println!("====================================================================");

        let ref_mol = mols_batch_next.remove(0);
        selected_mols.push(ref_mol.clone()); 
        let ref_eny = enys_batch_next.remove(0);
        selected_enys.push(ref_eny.clone()); 

        let mols_batch_curr = mols_batch_next.clone();
        mols_batch_next.clear();
        let enys_batch_curr = enys_batch_next.clone();
        enys_batch_next.clear();


        let (sender, receiver) = channel();

        let mols_size = mols_batch_curr.len();
        let thrds_num = rayon::current_num_threads();
        let mut  batch_size_thread = mols_size/thrds_num;
        let left_issue = mols_size%thrds_num;
        if left_issue > 0 {
            batch_size_thread += 1;
        }

        mols_batch_curr.par_chunks(batch_size_thread).zip(enys_batch_curr.par_chunks(batch_size_thread)).
            for_each_with(sender,|s, (sub_mols_batch, sub_enys_batch)| {

            let mut sub_mols_batch_next = vec![];
            let mut sub_enys_batch_next = vec![];
            sub_mols_batch.iter().zip(sub_enys_batch.iter()).for_each(|(mol,eny)| {
                //let (rot_mol, _, fmin_org, fmin) = similarity_optimization_lbfgs(mol, &ref_mol, displace, threshold);

                //if fmin > threshold {
                //    sub_mols_batch_next.push(rot_mol);
                //    sub_enys_batch_next.push(*eny);
                //}
                //if let Some(thread_index) = rayon::current_thread_index() {
                //    if thread_index == 0 {
                //        println!("thread 0 => fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);
                //    }
                //} else {
                //    println!("fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);

                //}
                let (flag, loss) = similarity_check(mol, &ref_mol, threshold, atm_type);
                if ! flag {
                    sub_mols_batch_next.push(mol.clone());
                    sub_enys_batch_next.push(*eny);
                };
                if let Some(thread_index) = rayon::current_thread_index() {
                    if thread_index == 0 {
                        println!("loss: {:16.8}", loss);
                    }
                } else {
                    println!("loss: {:16.8}", loss);
                }
            });
            s.send((sub_mols_batch_next, sub_enys_batch_next)).unwrap();
        });

        receiver.iter().for_each(|(sub_mols_batch_next, sub_enys_batch_next)| {
            mols_batch_next.extend(sub_mols_batch_next);
            enys_batch_next.extend(sub_enys_batch_next);
        });
    };

    if let Some(mol) = mols_batch_next.pop() {  
        selected_mols.push(mol);
    }
    if let Some(eny) = enys_batch_next.pop() {  
        selected_enys.push(eny);
    }

    (selected_mols, selected_enys)
}

pub fn filter_batch_serial(mols_batch: &[MatrixFull<f64>], _displace: f64, threshold: f64,atm_type: &[String]) -> Vec<MatrixFull<f64>> {
    //let mols_batch_curr = mols_batch.clone();

    let mut selected_mols: Vec<MatrixFull<f64>> = vec![];

    let mut mols_batch_next = mols_batch.to_vec();

    while mols_batch_next.len() > 1 {

        if let Some(thread_index) = rayon::current_thread_index() {
            let ref_mol = mols_batch_next.remove(0);
            selected_mols.push(ref_mol.clone()); 

            if thread_index == 0 {
                println!("====================================================================");
                println!(" Selected molecules: {:5}. Remaining molecules: {:5}", selected_mols.len(), mols_batch_next.len());
                println!("====================================================================");
            }

            let mols_batch_curr = mols_batch_next.clone();
            mols_batch_next.clear();

            mols_batch_curr.iter().for_each(|mol| {
                //let (rot_mot, _, fmin_org, fmin) = similarity_optimization_lbfgs(mol, &ref_mol, displace, threshold);

                //if fmin > threshold {
                //    mols_batch_next.push(rot_mot);
                //}
                //if thread_index == 0 {
                //    println!("fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);
                //}
                let (flag, loss) = similarity_check(mol, &ref_mol, threshold, atm_type);
                if ! flag {
                    mols_batch_next.push(mol.clone());
                    if thread_index == 0 {
                        println!("loss: {:16.8}", loss);
                    }
                }
            });
        } else {
            println!("====================================================================");
            println!(" Selected molecules: {:5}. Remaining molecules: {:5}", selected_mols.len(), mols_batch_next.len());
            println!("====================================================================");

            let ref_mol = mols_batch_next.remove(0);
            selected_mols.push(ref_mol.clone()); 

            let mols_batch_curr = mols_batch_next.clone();

            mols_batch_next.clear();

            mols_batch_curr.iter().for_each(|mol| {
                //let (rot_mot, _, fmin_org, fmin) = similarity_optimization_lbfgs(mol, &ref_mol, displace, threshold);
                //if fmin > threshold {
                //    mols_batch_next.push(rot_mot);
                //}
                //println!("fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);
                let (flag, loss) = similarity_check(mol, &ref_mol, threshold, atm_type);
                if ! flag {
                    mols_batch_next.push(mol.clone());
                    println!("loss: {:16.8}", loss);
                }
            });
        }
    }

    if let Some(mol) = mols_batch_next.pop() {  
        selected_mols.push(mol);
    }

    selected_mols

}

pub fn batch_check_rayon(mols_batch: &Vec<MatrixFull<f64>>, batch_size: usize, _displace: f64, threshold: f64, atm_type: &[String]) -> Vec<MatrixFull<f64>> {

    let mut mols_batch_next: Vec<MatrixFull<f64>> = vec![];

    let batch_size_curr = if batch_size > mols_batch.len() {mols_batch.len()} else {batch_size};

    let (sender, receiver) = channel();
    mols_batch.par_chunks(batch_size_curr).for_each_with(sender, |s, sub_mols_batch| {
        s.send(filter_batch_serial(sub_mols_batch, _displace, threshold, atm_type)).unwrap();
    });

    receiver.iter().for_each(|x| {
        mols_batch_next.extend(x);
    });

    mols_batch_next

}

pub fn batch_check(mols_batch: &Vec<MatrixFull<f64>>, energy_batch: &Option<Vec<f64>>, batch_size: usize, _displace: f64, threshold: f64, atm_type: &[String]) -> (Vec<MatrixFull<f64>>, Option<Vec<f64>>) {

    let mut mols_batch_next: Vec<MatrixFull<f64>> = vec![];

    let batch_size_curr = if batch_size > mols_batch.len() {mols_batch.len()} else {batch_size};
    let num_batchs = if mols_batch.len()%batch_size_curr == 0 {
        mols_batch.len()/batch_size_curr
    } else {
        mols_batch.len()/batch_size_curr + 1
    };

    match energy_batch {
        None => {

            mols_batch.chunks(batch_size_curr).enumerate().for_each(|(i,sub_mols_batch)| {
                println!("====================================================================");
                println!("Batch: {:5} of {:5}", i+1, num_batchs);
                println!("====================================================================");
                mols_batch_next.extend(filter_batch_rayon(sub_mols_batch, _displace, threshold, atm_type));
            });

            return (mols_batch_next, None)

        },
        Some(energy_batch) => {

            let mut enys_batch_next = vec![];

            mols_batch.chunks(batch_size_curr).zip(energy_batch.chunks(batch_size_curr)).enumerate()
            .for_each(|(i,(sub_mols_batch, sub_energy_batch))| {
                println!("====================================================================");
                println!("Batch: {:5} of {:5}", i+1, num_batchs);
                println!("====================================================================");
                let (selected_mols, selected_enys) = 
                    filter_batch_rayon_with_energy(sub_mols_batch, sub_energy_batch, _displace, threshold, atm_type);
                mols_batch_next.extend(selected_mols);
                enys_batch_next.extend(selected_enys);
            });

            return (mols_batch_next, Some(enys_batch_next))

        }

    }

}


pub fn format_xyz(mol: &MatrixFull<f64>, atm_type: &[String], energy: Option<f64>) -> anyhow::Result<String> {

    //let mut file = File::create(filename)?;
    let mut content = String::new();

    let atm_num = mol.size()[1];

    if atm_num != atm_type.len() {
        return Err(anyhow::anyhow!("Atom number mismatch"));
    }

    content.push_str(&format!("{}\n", atm_num));
    if let Some(ene) = energy {
        content.push_str(&format!(" energy={:-16.8} \n", ene))
    } else {
        content.push_str(&format!(" without energy \n"))
    };

    for i in 0..atm_num {
        content.push_str(&format!("{:3} {:16.8} {:16.8} {:16.8}\n", atm_type[i], mol[[0, i]], mol[[1, i]], mol[[2, i]]));
    }

    //file.write_all(content.as_bytes())?;

    Ok(content)
}
pub fn generate_output(out_file: &String, selected_mols: &Vec<MatrixFull<f64>>, selected_enys: &Option<Vec<f64>>, atm_type: &[String]) {
    let mut file = File::create(&out_file).unwrap();
    match selected_enys {
        None => {
            selected_mols.iter().for_each(|x| {
                let content = format_xyz(&x, &atm_type, None).unwrap();
                file.write(content.as_bytes()).unwrap();
            });
        },
        Some(selected_enys) => {
            selected_mols.iter().zip(selected_enys.iter()).for_each(|(x, e)| {
                let content = format_xyz(&x, &atm_type, Some(*e)).unwrap();
                file.write(content.as_bytes()).unwrap();
            });

        }
    }

}

pub fn filter_out_mode_0(
    mols_batch: &Vec<MatrixFull<f64>>, 
    energy_batch: &Option<Vec<f64>>, 
    atmtype_batch: &Vec<Vec<String>>,
    displace: f64,
    threshold: f64,
    batch_size: usize,
    fine_filter: bool) -> anyhow::Result<(Vec<MatrixFull<f64>>, Option<Vec<f64>>, Vec<String>)> {

    println!("Thresholds for similarity: {:16.8}", &threshold);
    if fine_filter {
        println!("Turn on the fine filter procedure");
    } else {
        println!("Turn off the fine filter procedure");

    }

    let atm_type = atmtype_batch[0].clone();

    let mut mols_batch_curr = mols_batch.clone();
    let mut enys_batch_curr = if let Some(enys_batch) = &energy_batch {
        if mols_batch.len() != enys_batch.len() {
            //return Err(anyhow::anyhow!("Geometry number mismatch"));
            println!("Energy number mismatch");
            None
        } else {
            energy_batch.clone()
        }
    } else {
        None
    };

    let mut mols_batch_next: Vec<MatrixFull<f64>> = vec![];
    let mut enys_batch_next = None;
    let mut iteration = true;
    let mut mols_size = mols_batch_curr.len();
    let mut i_iter = 0;

    let s_batch = init_timing();

    while iteration {
        let s_iter = init_timing();
        (mols_batch_next, enys_batch_next) = batch_check(&mols_batch_curr, &enys_batch_curr, batch_size, displace, threshold, &atm_type);

        let mols_size_curr = mols_batch_next.len();


        iteration = ! ((mols_size_curr == mols_size) || (mols_size-mols_size_curr < 10));
        if iteration {
            mols_size = mols_batch_next.len();
            mols_batch_curr = mols_batch_next.clone();
            enys_batch_curr = enys_batch_next.clone();
        }

        i_iter += 1;
        println!("====================================================================");
        println!(" After {:3} batch loops, filter out {:5} geometries", i_iter, mols_size_curr);
        timing(&s_iter, Some(&" This loop"));
        println!("====================================================================");
    }
    timing(&s_batch, Some(&"Block-based filter procedure"));

    let (selected_mols, selected_enys) = if ! fine_filter {
        (mols_batch_next, enys_batch_next)
    } else if batch_size < mols_batch_next.len() {
        let s_fine_filter = init_timing();
        let (filter_out_mols, filter_out_enys) = match enys_batch_next {
            None => {
                (filter_batch_rayon(&mols_batch_next, displace, threshold, &atm_type), None)
            },
            Some(enys_batch_next) => {
                let (tmp_mols, tmp_enys) = 
                    filter_batch_rayon_with_energy(&mols_batch_next, &enys_batch_next, displace, threshold, &atm_type);
                (tmp_mols, Some(tmp_enys))
            }
        };
        timing(&s_fine_filter, Some(&" Fine filter procedure"));
        (filter_out_mols, filter_out_enys)
    } else {
        (mols_batch_next, enys_batch_next)
    };

    println!("Selected molecules: {} -> {}", mols_batch.len(), selected_mols.len());
    Ok((selected_mols, selected_enys, atm_type))

}
pub fn filter_out_mode_1(mols_batch: &Vec<MatrixFull<f64>>, energy_batch: &Option<Vec<f64>>, atmtype_batch: &Vec<Vec<String>>, _displace: f64, threshold: f64, energy_threshold: f64) -> anyhow::Result<(Vec<MatrixFull<f64>>, Option<Vec<f64>>, Vec<String>)> {

    println!("Thresholds for similarity: {:16.8} and for energy: {:16.8}", &threshold, &energy_threshold);

    let atm_type = atmtype_batch[0].clone();

    let energy_batch = if let Some(e_list) = energy_batch {
        if e_list.len() == mols_batch.len() {
            e_list
        } else {
            panic!("Error: energy and geometry mismatch")
        }
    } else {
        panic!("Error: the energy information is needed for the filter mode 1")
    };


    let s_min = init_timing();
    let (min_mol, min_eny) = mols_batch.par_iter().zip(energy_batch.par_iter()).min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap()).unwrap();
    timing(&s_min, Some(&"Seeking for the minimal geometry"));

    let ref_mol = min_mol.clone();
    let ref_eny = min_eny.clone();
    println!("{}", format_xyz(&ref_mol, &atm_type, Some(ref_eny)).unwrap());

    let mut selected_batch: Vec<MatrixFull<f64>> = vec![ref_mol.clone()];
    let mut selected_enys: Vec<f64> = vec![*min_eny];

    //let selected_batch: Vec<(MatrixFull<f64>, f64)> = 
    let (sender, receiver) = channel();
    mols_batch.par_iter().zip(energy_batch.par_iter()).for_each_with(sender,|s, (x, e)| {
        //println!("{:16.8}, {:16.8}", e-ref_eny, e);
        if (e-ref_eny).abs() < energy_threshold {
            //let (rot_x, _, fmin_org, fmin) = similarity_optimization_lbfgs(x, &ref_mol, displace, threshold);
            //if fmin > threshold {
            //    println!("filter out a geometry with energy: {:16.8}, fmin_org: {:16.8}, fmin: {:16.8}", e, fmin_org, fmin);
            //    s.send((rot_x, *e)).unwrap();
            //} 
            let (flag, loss) = similarity_check(x, &ref_mol, threshold, &atm_type);
            if ! flag {
                println!("filter out a geometry with energy: {:16.8}, loss: {:16.8}", e, loss);
                s.send((x.clone(), *e)).unwrap();
            }
        } 
    });

    receiver.into_iter().for_each(|(rot_x,e)| {
        selected_batch.push(rot_x);
        selected_enys.push(e);
    });

    println!("Selected molecules: {} -> {}", mols_batch.len(), selected_batch.len()+1);

    Ok((selected_batch,Some(selected_enys), atm_type))
}

pub fn filter_out_mode_2(mols_batch: &Vec<MatrixFull<f64>>, energy_batch: &Option<Vec<f64>>, train_batch: &Vec<MatrixFull<f64>>, _displace: f64, threshold: f64, atm_type: &[String]) -> anyhow::Result<(Vec<MatrixFull<f64>>, Option<Vec<f64>>)> {

    //let num_trains = train_batch.len();

    println!("Thresholds for similarity: {:16.8}", &threshold);

    let s_rf = init_timing();

    let mut selected_mols: Vec<MatrixFull<f64>> = vec![];
    let mut selected_enys: Vec<f64> = vec![];


    match energy_batch {
        None => {
            let (sender, receiver) = channel();
            mols_batch.par_iter().enumerate().for_each_with(sender,|s, (i, x)| {
                let mut similar_flag = false;
                //let mut index = 0;
                train_batch.iter().try_for_each(|t| {
                    //let (_, _, _, fmin) = similarity_optimization_lbfgs(x, t, displace, threshold);
                    //if fmin <= threshold { 
                    //    similar_flag = true;
                    //    None
                    //} else {
                    //    Some(()) 
                    //}
                    let (flag, _) = similarity_check(x, t, threshold, atm_type);
                    if flag {
                        similar_flag = true;
                        None
                    } else {
                        Some(())
                    }
                });
                if ! similar_flag {
                    println!(" The {}th geometry is unique", i);
                    s.send(x).unwrap();
                }
            });
            receiver.into_iter().for_each(|x| {
                selected_mols.push((*x).clone());
            });
        },
        Some(energy_batch) => {
            let (sender, receiver) = channel();
            mols_batch.par_iter().zip(energy_batch.par_iter()).enumerate().for_each_with(sender,|s, (i,(x, e))| {
                let mut similar_flag = false;
                train_batch.iter().try_for_each(|t| {
                    //let (_, _, _, fmin) = similarity_optimization_lbfgs(x, t, displace, threshold);
                    //if fmin <= threshold { 
                    //    similar_flag = true;
                    //    None
                    //} else {
                    //    Some(()) 
                    //}
                    let (flag, _) = similarity_check(x, t, threshold, atm_type);
                    if flag {
                        similar_flag = true;
                        None
                    } else {
                        Some(())
                    }
                });

                if ! similar_flag {
                    println!(" The {}th geometry is unique", i);
                    s.send((x, e)).unwrap();
                }
            });
            receiver.into_iter().for_each(|(x,e)| {
                selected_mols.push((*x).clone());
                selected_enys.push(*e);
            });
        }

    }


    println!("Selected molecules: {} -> {}", mols_batch.len(), selected_mols.len());
    timing(&s_rf, Some(&"Second step in Mode 2"));

    match energy_batch {
        None => {
            Ok((selected_mols, None))
        },
        Some(_) => {
            Ok((selected_mols, Some(selected_enys)))
        }
    }

}


#[test]
fn test() {
    let ref_mol = MatrixFull::from_vec([3,6], 
        vec![-0.12595381,  0.95285406,  0.05294769, 
        0.        ,  0.        ,  0.        ,
        0.81185582, -0.13743935,  0.49853242,
        1.0359137 ,  0.24690515, -1.58758422,
        1.74925126,  0.6610954 , -2.09061483,
        1.6891999 ,  0.25419291, -2.95664034]
    ).unwrap();

    let mol = MatrixFull::from_vec([3,6],
        vec![
        -0.1145108 ,  0.93197491,  0.14567111 ,
         0.        ,  0.        ,  0.         ,
         0.80256949, -0.20205494,  0.50196594 ,
         0.97697016,  0.1908875 , -1.6331495  ,
         1.70793504,  0.6708945 , -2.06597271 ,
         1.66361823,  0.41985374, -2.97959794  
        ]
    ).unwrap();

    let psql = similarity_optimization_lbfgs(&mol, &ref_mol, 0.01, 0.1);

    println!("{:?}", &psql);
}

#[test]
fn regression() {
    let re_coord = Regex::new(
      r"(?P<elem>\w{1,2})\s+(?P<x>[\+-]?\d+.\d+)\s+(?P<y>[\+-]?\d+.\d+)\s+(?P<z>[\+-]?\d+.\d+)"
    ).unwrap();
    let re_energy = Regex::new(
      r"\s+energy=(?P<energy>[\+-]?\d+.\d+)\s+"
    ).unwrap();
    let re_pos_mom_for = Regex::new(
      r"(?P<elem>\w{1,2})\s+(?P<px>[\+-]?\d+.\d+)\s+(?P<py>[\+-]?\d+.\d+)\s+(?P<pz>[\+-]?\d+.\d+)\s+(?P<mx>[\+-]?\d+.\d+)\s+(?P<my>[\+-]?\d+.\d+)\s+(?P<mz>[\+-]?\d+.\d+)\s+(?P<fx>[\+-]?\d+.\d+)\s+(?P<fy>[\+-]?\d+.\d+)\s+(?P<fz>[\+-]?\d+.\d+)"
    ).unwrap();

    let match_line = "O       -0.83187344      -1.04911653      -1.21846370       0.85285893       0.29526215      -0.38053867       0.00653841       0.00055887       0.00042105".to_string();
    if re_coord.is_match(&match_line) {
        println!("coord line: {}", &match_line);
        let caps = re_coord.captures(&match_line).unwrap();
        let x = caps.name("x").unwrap().as_str().parse::<f64>().unwrap();
        let y = caps.name("y").unwrap().as_str().parse::<f64>().unwrap();
        let z = caps.name("z").unwrap().as_str().parse::<f64>().unwrap();
        let elem = caps.name("elem").unwrap().as_str();
        println!("{}: {:?}", &elem, &[x, y, z]);
    } else {
        println!("unknown line: {}", &match_line);
    }

    let match_line = "Properties=species:S:1:pos:R:3:momenta:R:3:forces:R:3 1S3W_xygjos_var13_mirror=T energy=-41.66744792223105 free_energy=-41.66744792223105 pbc=\"F F F\"".to_string();
    if re_energy.is_match(&match_line) {
        println!("energy line: {}", &match_line);
        let caps = re_energy.captures(&match_line).unwrap();
        let energy = caps.name("energy").unwrap().as_str().parse::<f64>().unwrap();
        println!("energy: {}", &energy);
    } else {
        println!("unknown line: {}", &match_line);
    };

    let match_line = "O    -0.83187344      -1.04911653      -1.21846370       0.85285893       0.29526215      -0.38053867       0.00653841       0.00055887       0.00042105".to_string();
    if re_pos_mom_for.is_match(&match_line) {
        println!("pos mom for line: {}", &match_line);
        let caps = re_pos_mom_for.captures(&match_line).unwrap();
        let px = caps.name("px").unwrap().as_str().parse::<f64>().unwrap();
        let py = caps.name("py").unwrap().as_str().parse::<f64>().unwrap();
        let pz = caps.name("pz").unwrap().as_str().parse::<f64>().unwrap();
        let mx = caps.name("mx").unwrap().as_str().parse::<f64>().unwrap();
        let my = caps.name("my").unwrap().as_str().parse::<f64>().unwrap();
        let mz = caps.name("mz").unwrap().as_str().parse::<f64>().unwrap();
        let fx = caps.name("fx").unwrap().as_str().parse::<f64>().unwrap();
        let fy = caps.name("fy").unwrap().as_str().parse::<f64>().unwrap();
        let fz = caps.name("fz").unwrap().as_str().parse::<f64>().unwrap();
        let elem = caps.name("elem").unwrap().as_str();
        println!("{}: {:?}", &elem, &[px, py, pz, mx, my, mz, fx, fy, fz]);
    } else {
        println!("unknown line: {}", &match_line);
    }

}

#[test]
fn test_fm() {
    println!("{}", formated_element_name(&String::from("he")));
    let chunk_size = 5_usize;
    let mols_batch = (0..24).map(|x| x as f64).collect::<Vec<f64>>();
    mols_batch.par_chunks(chunk_size).for_each(|x| {
        println!("thread {}: {:?}", rayon::current_thread_index().unwrap(), x);
    })
}


pub fn filter_out_mode_0_with_index(
    mols_batch: &Vec<MatrixFull<f64>>, 
    atmtype_batch: &Vec<Vec<String>>,
    displace: f64,
    threshold: f64,
    batch_size: usize,
    fine_filter: bool) -> anyhow::Result<(Vec<MatrixFull<f64>>, Vec<usize>, Vec<String>)> {

    println!("Thresholds for similarity: {:16.8}", &threshold);
    if fine_filter {
        println!("Turn on the fine filter procedure");
    } else {
        println!("Turn off the fine filter procedure");

    }

    let atm_type = atmtype_batch[0].clone();

    let mut mols_batch_curr = mols_batch.clone();
    let mut indx_batch_curr = (0..mols_batch_curr.len()).collect::<Vec<usize>>();

    let mut mols_batch_next: Vec<MatrixFull<f64>> = vec![];
    let mut indx_batch_next: Vec<usize> = vec![];

    let mut iteration = true;
    let mut mols_size = mols_batch_curr.len();
    let mut i_iter = 0;

    let s_batch = init_timing();

    while iteration {
        let s_iter = init_timing();
        (mols_batch_next, indx_batch_next) = batch_check_with_index(&mols_batch_curr, &indx_batch_curr, batch_size, displace, threshold, &atm_type);

        let mols_size_curr = mols_batch_next.len();


        iteration = ! ((mols_size_curr == mols_size) || (mols_size-mols_size_curr < 10));
        if iteration {
            mols_size = mols_batch_next.len();
            mols_batch_curr = mols_batch_next.clone();
            indx_batch_curr = indx_batch_next.clone();
        }

        i_iter += 1;
        println!("====================================================================");
        println!(" After {:3} batch loops, filter out {:5} geometries", i_iter, mols_size_curr);
        timing(&s_iter, Some(&" This loop"));
        println!("====================================================================");
    }
    timing(&s_batch, Some(&"Block-based filter procedure"));

    let (selected_mols, selected_indx) = if ! fine_filter {
        (mols_batch_next, indx_batch_next)
    } else if batch_size < mols_batch_next.len() {
        let s_fine_filter = init_timing();

        let (filter_out_mols, filter_out_indx) = filter_batch_rayon_with_index(&mols_batch_next, &indx_batch_next, displace, threshold, &atm_type);

        timing(&s_fine_filter, Some(&" Fine filter procedure"));
        (filter_out_mols, filter_out_indx)
    } else {
        (mols_batch_next, indx_batch_next)
    };

    println!("Selected molecules: {} -> {}", mols_batch.len(), selected_mols.len());
    Ok((selected_mols, selected_indx, atm_type))

}

pub fn batch_check_with_index(mols_batch: &Vec<MatrixFull<f64>>, index_batch: &Vec<usize>, batch_size: usize, displace: f64, threshold: f64, atm_type: &[String]) -> (Vec<MatrixFull<f64>>, Vec<usize>) {

    let mut mols_batch_next: Vec<MatrixFull<f64>> = vec![];
    let mut indx_batch_next: Vec<usize> = vec![];

    let batch_size_curr = if batch_size > mols_batch.len() {mols_batch.len()} else {batch_size};
    let num_batchs = if mols_batch.len()%batch_size_curr == 0 {
        mols_batch.len()/batch_size_curr
    } else {
        mols_batch.len()/batch_size_curr + 1
    };

    mols_batch.chunks(batch_size_curr).zip(index_batch.chunks(batch_size_curr)).enumerate()
    .for_each(|(i,(sub_mols_batch, sub_indx_batch))| {
        println!("====================================================================");
        println!("Batch: {:5} of {:5}", i+1, num_batchs);
        println!("====================================================================");
        let (selected_mols, selected_indx) = filter_batch_rayon_with_index(sub_mols_batch, sub_indx_batch, displace, threshold, atm_type);
        mols_batch_next.extend(selected_mols);
        indx_batch_next.extend(selected_indx);
    });

    return (mols_batch_next, indx_batch_next)

}

pub fn filter_batch_rayon_with_index(mols_batch: &[MatrixFull<f64>], index_batch: &[usize], _displace: f64, threshold: f64, atm_type: &[String]) -> (Vec<MatrixFull<f64>>, Vec<usize>) {

    //let thread_number = rayon::current_num_threads();

    let mut selected_mols: Vec<MatrixFull<f64>> = vec![];
    let mut selected_indx: Vec<usize> = vec![];

    let mut mols_batch_next = mols_batch.to_vec();
    let mut indx_batch_next = index_batch.to_vec();

    while mols_batch_next.len() > 1 {

        println!("====================================================================");
        println!(" Selected molecules: {:5}. Remaining molecules: {:5}", selected_mols.len(), mols_batch_next.len());
        println!("====================================================================");

        let ref_mol = mols_batch_next.remove(0);
        selected_mols.push(ref_mol.clone()); 
        let ref_idx = indx_batch_next.remove(0);
        selected_indx.push(ref_idx.clone()); 

        let mols_batch_curr = mols_batch_next.clone();
        mols_batch_next.clear();
        let indx_batch_curr = indx_batch_next.clone();
        indx_batch_next.clear();


        let (sender, receiver) = channel();

        let mols_size = mols_batch_curr.len();
        let thrds_num = rayon::current_num_threads();
        let mut  batch_size_thread = mols_size/thrds_num;
        let left_issue = mols_size%thrds_num;
        if left_issue > 0 {
            batch_size_thread += 1;
        }

        mols_batch_curr.par_chunks(batch_size_thread).zip(indx_batch_curr.par_chunks(batch_size_thread)).
            for_each_with(sender,|s, (sub_mols_batch, sub_indx_batch)| {

            let mut sub_mols_batch_next = vec![];
            let mut sub_indx_batch_next = vec![];
            sub_mols_batch.iter().zip(sub_indx_batch.iter()).for_each(|(mol,idx)| {
                //let (rot_mol, _, fmin_org, fmin) = similarity_optimization_lbfgs(mol, &ref_mol, displace, threshold);
                //if fmin > threshold {
                //    sub_mols_batch_next.push(rot_mol);
                //    sub_indx_batch_next.push(*idx);
                //}
                //if let Some(thread_index) = rayon::current_thread_index() {
                //    if thread_index == 0 {
                //        println!("thread 0 => fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);
                //    }
                //} else {
                //    println!("fmin_org: {:16.8}, fmin: {:16.8}", fmin_org, fmin);

                //}
                let (flag, loss) = similarity_check(mol, &ref_mol, threshold, &atm_type);
                if ! flag {
                    sub_mols_batch_next.push(mol.clone());
                    sub_indx_batch_next.push(*idx);
                }
                if let Some(thread_index) = rayon::current_thread_index() {
                    if thread_index == 0 {
                        println!("thread 0 => loss: {:16.8}", loss);
                    }
                } else {
                    println!("thread 0 => loss: {:16.8}", loss);
                }
            });
            s.send((sub_mols_batch_next, sub_indx_batch_next)).unwrap();
        });

        receiver.iter().for_each(|(sub_mols_batch_next, sub_indx_batch_next)| {
            mols_batch_next.extend(sub_mols_batch_next);
            indx_batch_next.extend(sub_indx_batch_next);
        });
    };

    if let Some(mol) = mols_batch_next.pop() {  
        selected_mols.push(mol);
    }
    if let Some(idx) = indx_batch_next.pop() {  
        selected_indx.push(idx);
    }

    (selected_mols, selected_indx)
}

pub fn filter_out_mode_1_with_index(mols_batch: &Vec<MatrixFull<f64>>, energy_batch: &Option<Vec<f64>>, atmtype_batch: &Vec<Vec<String>>, _displace: f64, threshold: f64, energy_threshold: f64) -> anyhow::Result<(Vec<MatrixFull<f64>>, Vec<usize>, Vec<String>)> {

    println!("Thresholds for similarity: {:16.8} and for energy: {:16.8}", &threshold, &energy_threshold);

    let atm_type = atmtype_batch[0].clone();

    let energy_batch = if let Some(e_list) = energy_batch {
        if e_list.len() == mols_batch.len() {
            e_list
        } else {
            panic!("Error: energy and geometry mismatch")
        }
    } else {
        panic!("Error: the energy information is needed for the filter mode 1")
    };

    let s_min = init_timing();
    let ((min_idx, min_mol), min_eny) = mols_batch.par_iter().enumerate().zip(energy_batch.par_iter())
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap()).unwrap();
    timing(&s_min, Some(&"Seeking for the minimal geometry"));

    let ref_mol = min_mol.clone();
    let ref_eny = min_eny.clone();
    println!("{}", format_xyz(&ref_mol, &atm_type, Some(ref_eny)).unwrap());

    let mut selected_batch: Vec<MatrixFull<f64>> = vec![ref_mol.clone()];
    //let mut selected_enys: Vec<f64> = vec![];
    let mut selected_index: Vec<usize> = vec![min_idx];

    //let selected_batch: Vec<(MatrixFull<f64>, f64)> = 
    let (sender, receiver) = channel();
    mols_batch.par_iter().zip(energy_batch.par_iter()).enumerate().for_each_with(sender,|s, (i, (x, e))| {
        //println!("{:16.8}, {:16.8}", e-ref_eny, e);
        if (e-ref_eny).abs() < energy_threshold {
            //let (rot_x, _, fmin_org, fmin) = similarity_optimization_lbfgs(x, &ref_mol, displace, threshold);
            //if fmin > threshold {
            //    println!("filter out a geometry with energy: {:16.8}, fmin_org: {:16.8}, fmin: {:16.8}", e, fmin_org, fmin);
            //    //s.send((rot_x, *e, i)).unwrap();
            //    s.send((rot_x, i)).unwrap();
            //} 
            let (flag, loss) = similarity_check(x, &ref_mol, threshold, &atm_type);
            if ! flag {
                println!("filter out a geometry with energy: {:16.8}, loss: {:16.8}", e, loss);
                s.send((x.clone(), i)).unwrap();
            }
        } 
    });

    //receiver.into_iter().for_each(|(rot_x,e, i)| {
    receiver.into_iter().for_each(|(rot_x,i)| {
        selected_batch.push(rot_x);
        //selected_enys.push(e);
        selected_index.push(i);
    });

    println!("Selected molecules: {} -> {}", mols_batch.len(), selected_batch.len()+1);

    Ok((selected_batch,selected_index, atm_type))
}

pub fn filter_out_mode_2_with_index(mols_batch: &Vec<MatrixFull<f64>>, indx_batch: &Vec<usize>, train_batch: &Vec<MatrixFull<f64>>, _displace: f64, threshold: f64, atm_type: &[String]) -> anyhow::Result<(Vec<MatrixFull<f64>>, Vec<usize>)> {

    //let num_trains = train_batch.len();

    println!("Thresholds for similarity: {:16.8}", &threshold);

    let s_rf = init_timing();

    let mut selected_mols: Vec<MatrixFull<f64>> = vec![];
    let mut selected_indx: Vec<usize> = vec![];


    let (sender, receiver) = channel();
    mols_batch.par_iter().zip(indx_batch.par_iter()).for_each_with(sender,|s, (x, i)| {
        let mut similar_flag = false;
        train_batch.iter().try_for_each(|t| {
            //let (_, _, _, fmin) = similarity_optimization_lbfgs(x, t, displace, threshold);
            //if fmin <= threshold { 
            //    similar_flag = true;
            //    None
            //} else {
            //    Some(()) 
            //}
            let (flag, _) = similarity_check(x, t, threshold, atm_type);
            if flag {
                similar_flag = true;
                None
            } else {
                Some(())
            }
        });

        if ! similar_flag {
            println!(" The {}th geometry is unique", i);
            s.send((x, i)).unwrap();
        }
    });
    receiver.into_iter().for_each(|(x,i)| {
        selected_mols.push((*x).clone());
        selected_indx.push(*i);
    });



    println!("Selected molecules: {} -> {}", mols_batch.len(), selected_mols.len());
    timing(&s_rf, Some(&"Second step in Mode 2"));

    Ok((selected_mols,selected_indx))

}

pub fn generate_output_with_index(out_file: &String, selected_mols: &Vec<MatrixFull<f64>>, selected_indx: &Vec<usize>, atm_type: &[String], energy_batch: &Option<Vec<f64>>, momfor_batch: &Option<Vec<MatrixFull<f64>>>) {
    let mut file = File::create(&out_file).unwrap();

    match (energy_batch, momfor_batch) {
        (None, None) => {
            selected_mols.iter().for_each(|x| {
                let content = format_xyz_v02(&x, &atm_type, None, None).unwrap();
                file.write(content.as_bytes()).unwrap();
            });
        },
        (Some(enys), Some(mafs)) => {
            selected_mols.iter().zip(selected_indx.iter()).for_each(|(x, i)| {
                let e = enys.get(*i).unwrap();
                let mf = mafs.get(*i).unwrap();
                let content = format_xyz_v02(&x, &atm_type, Some(*e), Some(mf)).unwrap();
                file.write(content.as_bytes()).unwrap();
            });
        },
        (None, Some(mafs)) => {
            selected_mols.iter().zip(selected_indx.iter()).for_each(|(x, i)| {
                //let e = enys.get(*i).unwrap();
                let mf = mafs.get(*i).unwrap();
                let content = format_xyz_v02(&x, &atm_type, None, Some(mf)).unwrap();
                file.write(content.as_bytes()).unwrap();
            });
        },
        (Some(enys), None) => {
            selected_mols.iter().zip(selected_indx.iter()).for_each(|(x, i)| {
                let e = enys.get(*i).unwrap();
                //let mf = mafs.get(*i).unwrap();
                let content = format_xyz_v02(&x, &atm_type, Some(*e), None).unwrap();
                file.write(content.as_bytes()).unwrap();
            });
        },
    }

}

pub fn format_xyz_v02(mol: &MatrixFull<f64>, atm_type: &[String], energy: Option<f64>, mf: Option<&MatrixFull<f64>>) -> anyhow::Result<String> {

    let mut content = String::new();

    let atm_num = mol.size()[1];

    if atm_num != atm_type.len() {
        return Err(anyhow::anyhow!("Atom number mismatch"));
    }

    content.push_str(&format!("{}\n", atm_num));
    if let Some(ene) = energy {
        content.push_str(&format!(" energy={:-16.8} \n", ene))
    } else {
        content.push_str(&format!(" without energy \n"))
    };

    if let Some(mf) = mf {
        for i in 0..atm_num {
            content.push_str(&format!("{:3} {:16.8} {:16.8} {:16.8} {:16.8} {:16.8} {:16.8} {:16.8} {:16.8} {:16.8}\n", 
                atm_type[i], mol[[0, i]], mol[[1, i]], mol[[2, i]], mf[[0, i]], mf[[1, i]], mf[[2, i]], mf[[3, i]], mf[[4, i]], mf[[5, i]]));
        }
    } else {
        for i in 0..atm_num {
            content.push_str(&format!("{:3} {:16.8} {:16.8} {:16.8}\n", atm_type[i], mol[[0, i]], mol[[1, i]], mol[[2, i]]));
        }
    }

    Ok(content)
}

fn main() -> anyhow::Result<()>{

    let s_time = init_timing();
    let mode: usize = parse_input().get_one::<String>("mode").unwrap_or(&String::from("0")).parse().unwrap();
    println!("Action mode: {}", mode);

    let ctrl_file = parse_input().get_one::<String>("input_file").unwrap_or(&String::from("geometries.xyz")).to_string();
    let out_file = parse_input().get_one::<String>("output_file").unwrap_or(&String::from("filter_out.xyz")).to_string();
    println!("Input geometries in {} and output geometries in {}", &ctrl_file, &out_file);

    let displace = parse_input().get_one::<String>("displace").unwrap_or(&String::from("0.001")).parse().unwrap();

    let num_threads: usize = parse_input().get_one::<String>("num_threads").unwrap_or(&String::from("8")).parse().unwrap();
    rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global().unwrap();
    println!("Rayon parallelism with {} threads", rayon::current_num_threads());

    let s_rf = init_timing();
    let (mut mols_batch, energy_batch, atmtype_batch, momfor_batch) = import_xyz_in_batch(&ctrl_file).unwrap();
    println!("mols_batch: ({},{})", &ctrl_file, mols_batch.len());
    timing(&s_rf, Some(&"Importing geometries"));

    let atm_type = atmtype_batch[0].clone();
    let mass_list = get_mass(&atm_type);
    mols_batch.iter_mut().for_each(|x| {
        //shift_wrt_ref_atom(x, 1)
        shift_to_mass_center(x, &mass_list)
    });

    let mut selected_mols = vec![];
    //let mut selected_enys= None;
    let mut selected_indx= vec![];
    //let mut atm_type = vec![];
    if mode == 0 {

        let threshold = parse_input().get_one::<String>("similarity_threshold").unwrap_or(&String::from("0.2")).parse().unwrap();
        let batch_size: usize = parse_input().get_one::<String>("block_size").unwrap_or(&String::from("100")).parse().unwrap();
        let fine_filter = parse_input().get_flag("fine_filter");

        (selected_mols, selected_indx, _) = filter_out_mode_0_with_index(&mols_batch, &atmtype_batch, displace, threshold, batch_size, fine_filter).unwrap();

    } else if mode == 1 {
        let threshold = parse_input().get_one::<String>("similarity_threshold").unwrap_or(&String::from("0.2")).parse().unwrap();
        let energy_threshold: f64 = parse_input().get_one::<String>("energy_threshold").unwrap_or(&String::from("0.001")).parse().unwrap();

        (selected_mols, selected_indx, _) = filter_out_mode_1_with_index(&mols_batch, &energy_batch, &atmtype_batch, displace, threshold, energy_threshold).unwrap();

    } else if mode == 2 {
        let threshold = parse_input().get_one::<String>("similarity_threshold").unwrap_or(&String::from("0.2")).parse().unwrap();
        let batch_size: usize = parse_input().get_one::<String>("block_size").unwrap_or(&String::from("100")).parse().unwrap();
        let fine_filter = parse_input().get_flag("fine_filter");

        (selected_mols, selected_indx, _) = filter_out_mode_0_with_index(&mols_batch,  &atmtype_batch, displace, threshold, batch_size, fine_filter).unwrap();

        let train_file = parse_input().get_one::<String>("training_file").unwrap_or(&String::from("training.xyz")).to_string(); 
        let (mut train_mols_batch, _, _, _ )= import_xyz_in_batch(&train_file).unwrap_or((vec![], None, vec![], None));
        println!("Trainning set: ({}, {})", &train_file, &train_mols_batch.len());
        train_mols_batch.iter_mut().for_each(|x| {
            //shift_wrt_ref_atom(x, 1)
            shift_to_mass_center(x, &mass_list)
        });
        (selected_mols, selected_indx) = filter_out_mode_2_with_index(&selected_mols, &selected_indx, &train_mols_batch, displace, threshold, &atm_type).unwrap();

    } else if mode == 3 {

        let number_threshold = parse_input().get_one::<String>("number_threshold").unwrap_or(&String::from("1000")).parse().unwrap();
        let batch_size: usize = parse_input().get_one::<String>("block_size").unwrap_or(&String::from("100")).parse().unwrap();
        let fine_filter = parse_input().get_flag("fine_filter");

        let train_file = parse_input().get_one::<String>("training_file").unwrap_or(&String::from("training.xyz")).to_string(); 
        let (mut train_mols_batch, _ , _, _) = import_xyz_in_batch(&train_file).unwrap_or((vec![], None, vec![], None));
        println!("Trainning set: ({}, {})", &train_file, &train_mols_batch.len());
        train_mols_batch.iter_mut().for_each(|x| {
            //shift_wrt_ref_atom(x, 1)
            shift_to_mass_center(x, &mass_list)
        });

        let mut threshold = 1.0;

        while (selected_mols.len() < number_threshold) && selected_mols.len() < mols_batch.len() {

            //let (mut tmp_mols, mut tmp_enys, _) = filter_out_mode_0(&mols_batch, &energy_batch, &atmtype_batch, displace, threshold, batch_size, fine_filter).unwrap();
            //(tmp_mols, tmp_enys) = filter_out_mode_2(&tmp_mols, &tmp_enys, &train_mols_batch, displace, threshold).unwrap();
            //selected_mols.append(&mut tmp_mols);
            //if let Some(mut tmp_enys) = tmp_enys {
            //    if let Some(selected_enys) = &mut selected_enys {
            //        selected_enys.append(&mut tmp_enys);
            //    }
            //}
            //train_mols_batch=selected_mols.clone();
            (selected_mols, selected_indx, _) = filter_out_mode_0_with_index(&mols_batch,  &atmtype_batch, displace, threshold, batch_size, fine_filter).unwrap();
            

            threshold *= 0.5;

            println!("====================================================================");
            println!(" For mode 3, filter out {:5} geometries in total after this round", selected_mols.len());
            println!("====================================================================");
        }

        if selected_mols.len() > number_threshold {
            selected_mols.truncate(number_threshold);
            selected_indx.truncate(number_threshold);
        }

    }

    generate_output_with_index(&out_file, &selected_mols, &selected_indx, &atm_type, &energy_batch, &momfor_batch);
    timing(&s_time, Some(&"Total time consumption"));

    Ok(())

}