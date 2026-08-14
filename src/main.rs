use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let arg_strs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
    if let Err(e) = treecode::treecode::run(&arg_strs) {
        match e {
            treecode::error::TreeError::Help => std::process::exit(0),
            _ => {
                eprintln!("{}", e);
                std::process::exit(1);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn binary_entry_runs_simulation() {
        let args = ["treecode", "nbody=30", "tstop=0.01", "dtout=0.005"];
        let tree = treecode::treecode::run(&args).unwrap();
        assert!(tree.nstep > 0, "simulation should advance steps");
    }
}
