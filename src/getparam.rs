use crate::error::{Result, TreeError};

const DEFPARAM: i32 = 0o1;
const REQPARAM: i32 = 0o2;
const ARGPARAM: i32 = 0o4;

#[derive(Clone, Debug)]
struct Param {
    name: String,
    value: String,
    comment: Option<String>,
    flags: i32,
}

/// Owns the parsed command-line parameters. This replaces the C `static`
/// parameter table and the `PROGNAME` global: a [`Config`] is built once by
/// [`Config::initparam`] (the C `initparam`, now an associated function) and
/// then queried through the typed getters that mirror the original C API
/// (`getparam`, `getparamstat`, `getiparam`, `getdparam`, `getbparam`).
#[derive(Debug)]
pub struct Config {
    progname: String,
    params: Vec<Param>,
}

impl Default for Config {
    fn default() -> Self {
        Config::new()
    }
}

impl Config {
    pub fn new() -> Self {
        Config {
            progname: String::new(),
            params: Vec::new(),
        }
    }

    pub fn initparam(argv: &[&str], defv: &[&str]) -> Result<Config> {
        let mut config = Config::new();
        config.progname = argv[0].to_string();

        config.params.push(Param {
            name: "argv0".to_string(),
            value: argv[0].to_string(),
            comment: None,
            flags: ARGPARAM,
        });

        let mut current_comment: Option<String> = None;
        defv.iter().for_each(|&entry| {
            if entry.starts_with(';') {
                current_comment = Some(entry.strip_prefix(';').unwrap_or(entry).to_string());
            } else {
                let (name, value) = parse_name_value(entry);
                let flags = if value == "???" {
                    DEFPARAM | REQPARAM
                } else {
                    DEFPARAM
                };
                let comment = current_comment.take();
                config.params.push(Param {
                    name,
                    value,
                    comment,
                    flags,
                });
            }
        });

        if argv.len() > 1 && (argv[1] == "-clue" || argv[1] == "-help") {
            let p = &config.params[0];
            if argv[1] == "-clue" {
                print!("{}", p.value);
                config.params[1..].iter().for_each(|pp| {
                    print!(" {}={}", pp.name, pp.value);
                });
                println!();
            } else {
                println!("{}", p.value);
                config.params[1..].iter().for_each(|pp| {
                    let item = format!("  {}={}", pp.name, pp.value);
                    if let Some(ref c) = pp.comment {
                        if item.len() < 32 {
                            println!("{:<32}  {}", item, c);
                        } else {
                            println!("{}\t\t\t\t  {}", item, c);
                        }
                    } else {
                        println!("{}", item);
                    }
                });
            }
            return Err(TreeError::Help);
        }

        let mut scanpos = true;
        let mut pidx = 0;
        argv[1..].iter().try_for_each(|&arg| {
            if let Some((name, value)) = parse_name_value_opt(arg) {
                scanpos = false;
                if let Some(pp) = config.params.iter_mut().find(|p| p.name == name) {
                    if pp.flags & ARGPARAM != 0 {
                        return Err(TreeError::ParamDuplicated {
                            prog: argv[0].to_string(),
                            name,
                        });
                    }
                    pp.value = value;
                    pp.flags = (pp.flags & !DEFPARAM) | ARGPARAM;
                } else {
                    return Err(TreeError::UnknownParam {
                        prog: argv[0].to_string(),
                        name,
                    });
                }
            } else if scanpos {
                pidx += 1;
                if pidx >= config.params.len() {
                    return Err(TreeError::TooManyArgs(argv[0].to_string()));
                }
                config.params[pidx].value = arg.to_string();
                config.params[pidx].flags = (config.params[pidx].flags & !DEFPARAM) | ARGPARAM;
            } else {
                return Err(TreeError::NamelessArg {
                    prog: argv[0].to_string(),
                    arg: arg.to_string(),
                });
            }
            Ok(())
        })?;

        let mut needarg = false;
        config.params[1..].iter().for_each(|pp| {
            if (pp.flags & REQPARAM != 0) && (pp.flags & DEFPARAM != 0) {
                needarg = true;
            }
        });
        if needarg {
            eprint!("Usage: {}", config.params[0].value);
            config.params[1..].iter().for_each(|pp| {
                if pp.flags & REQPARAM != 0 {
                    eprint!(" {}=???", pp.name);
                }
            });
            eprintln!(": required arguments missing");
            return Err(TreeError::MissingRequiredParam);
        }

        Ok(config)
    }

    pub fn getparam(&self, name: &str) -> Result<String> {
        if let Some(pp) = self.params.iter().find(|p| p.name == name) {
            return Ok(pp.value.clone());
        }
        if name == "argv0" {
            return Ok(self.progname.clone());
        }
        Err(TreeError::ParamNotAvailable(name.to_string()))
    }

    pub fn getparamstat(&self, name: &str) -> i32 {
        if let Some(pp) = self.params.iter().find(|p| p.name == name) {
            return pp.flags;
        }
        0
    }

    pub fn getiparam(&self, name: &str) -> Result<i32> {
        let val_str = self.getparam(name)?;
        let val: i64 = val_str.parse().unwrap_or_else(|_| {
            let trimmed = val_str.trim_end_matches(|c: char| c.is_alphabetic());
            let suffix = val_str[trimmed.len()..].to_string();
            let base: i64 = trimmed.parse().unwrap_or(0);
            match suffix.as_str() {
                "k" | "K" => base * 1024,
                "m" | "M" => base * 1024 * 1024,
                _ => base,
            }
        });
        Ok(val as i32)
    }

    pub fn getdparam(&self, name: &str) -> Result<f64> {
        let val_str = self.getparam(name)?;
        if let Some((n, d)) = val_str.split_once('/') {
            let n: f64 = n.parse().unwrap_or(0.0);
            let d: f64 = d.parse().unwrap_or(1.0);
            Ok(n / d)
        } else {
            val_str
                .parse()
                .map_err(|_| TreeError::ParamNotAvailable(name.to_string()))
        }
    }

    pub fn getbparam(&self, name: &str) -> Result<bool> {
        let val = self.getparam(name)?;
        let first = val.chars().next().unwrap_or(' ');
        match first {
            't' | 'T' | 'y' | 'Y' | '1' => Ok(true),
            'f' | 'F' | 'n' | 'N' | '0' => Ok(false),
            _ => Err(TreeError::BadBoolParam {
                name: name.to_string(),
                value: val,
            }),
        }
    }
}

fn parse_name_value(s: &str) -> (String, String) {
    if let Some(eq) = s.find('=') {
        (s[..eq].to_string(), s[eq + 1..].to_string())
    } else {
        (s.to_string(), String::new())
    }
}

fn parse_name_value_opt(s: &str) -> Option<(String, String)> {
    let s = s.strip_prefix(['<', '>']).unwrap_or(s);
    s.find('=')
        .map(|eq| (s[..eq].to_string(), s[eq + 1..].to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn getparam_after_init() -> Result<()> {
        let cfg = Config::initparam(&["prog", "x=5"], &["x=0"])?;
        assert_eq!(cfg.getparam("x")?, "5");
        assert_eq!(cfg.getparamstat("x") & ARGPARAM, ARGPARAM);
        Ok(())
    }

    #[test]
    fn getparam_argv0() -> Result<()> {
        let cfg = Config::initparam(&["hello"], &[])?;
        assert_eq!(cfg.getparam("argv0")?, "hello");
        Ok(())
    }

    #[test]
    fn parse_helpers() {
        assert_eq!(
            parse_name_value("name=value"),
            ("name".to_string(), "value".to_string())
        );
        assert_eq!(
            parse_name_value("bare"),
            ("bare".to_string(), String::new())
        );
        assert_eq!(
            parse_name_value_opt("<x=1"),
            Some(("x".to_string(), "1".to_string()))
        );
        assert_eq!(parse_name_value_opt("noeq"), None);
    }
}
