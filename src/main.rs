#![forbid(clippy::unwrap_used)]
#![forbid(clippy::expect_used)]

use std::{env, process::exit};

fn main() {
    let params: Vec<String> = env::args().skip(1).collect();
    let result = match treecode::Simulation::new(params) {
        Ok(mut sim) => sim.run(),
        Err(e) => Err(e),
    };
    if let Err(e) = result {
        match e {
            treecode::error::TreeError::Help => exit(0),
            _ => {
                eprintln!("{}", e);
                exit(1);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn binary_entry_runs_simulation() -> Result<(), treecode::error::TreeError> {
        let args = ["treecode", "nbody=30", "tstop=0.01", "dtout=0.005"];
        let tree = treecode::treecode::run(&args)?;
        assert!(tree.nstep > 0, "simulation should advance steps");
        Ok(())
    }
}
