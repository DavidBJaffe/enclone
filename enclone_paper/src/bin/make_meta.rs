// Make a meta file for the test data.

use enclone_base::expand_integer_ranges;
use enclone_core::test_def::*;

fn main() {
    println!("bcr,donor,origin");
    let test1 = TEST1.split(':').collect::<Vec<&str>>();
    for (i, t) in test1.iter().enumerate() {
        let x = expand_integer_ranges(*t);
        let x = x.split(',').collect::<Vec<&str>>();
        for id in x.iter() {
            println!("{id},d1,o{}", i + 1);
        }
    }
    let test2 = TEST2.split(':').collect::<Vec<&str>>();
    for (i, t) in test2.iter().enumerate() {
        let x = expand_integer_ranges(*t);
        let x = x.split(',').collect::<Vec<&str>>();
        for id in x.iter() {
            println!("{id},d2,o{}", i + 1);
        }
    }
    let test3 = TEST3.split(':').collect::<Vec<&str>>();
    for (i, t) in test3.iter().enumerate() {
        let x = expand_integer_ranges(*t);
        let x = x.split(',').collect::<Vec<&str>>();
        for id in x.iter() {
            println!("{id},d3,o{}", i + 1);
        }
    }
    let test4 = TEST4.split(':').collect::<Vec<&str>>();
    for (i, t) in test4.iter().enumerate() {
        let x = expand_integer_ranges(*t);
        let x = x.split(',').collect::<Vec<&str>>();
        for id in x.iter() {
            println!("{id},d4,o{}", i + 1);
        }
    }
}
