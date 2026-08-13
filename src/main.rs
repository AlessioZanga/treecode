use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let arg_strs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
    treecode::treecode::run(&arg_strs);
}

#[cfg(test)]
mod tests {
    #[test]
    fn binary_entry_runs_simulation() {
        let args = ["treecode", "nbody=30", "tstop=0.01", "dtout=0.005"];
        treecode::treecode::run(&args);
        let nstep = unsafe { treecode::types::nstep };
        assert!(nstep > 0, "simulation should advance steps");
    }
}
