use std::io::Write;

pub type Result<T> = std::result::Result<T, TreeError>;

#[derive(Debug, thiserror::Error)]
pub enum TreeError {
    #[error("absurd value for nbody: {0}")]
    AbsurdNbody(i32),

    #[error("inputdata: ndim = {got}; expected {expected}")]
    BadNdim { got: i32, expected: usize },

    #[error("incompatible options bh86 and sw94")]
    IncompatibleOptions,

    #[error("hackcofm: tree structure error")]
    TreeStructure,

    #[error("walktree: active list overflow")]
    ActiveListOverflow,

    #[error("loadbody: two bodies have same position")]
    CoincidentBodies,

    #[error("walktree: recursion terminated with cell")]
    RecursionTerminated,

    #[error("{prog}: parameter {name} duplicated")]
    ParamDuplicated { prog: String, name: String },

    #[error("{prog}: parameter {name} unknown")]
    UnknownParam { prog: String, name: String },

    #[error("{0}: too many arguments")]
    TooManyArgs(String),

    #[error("{prog}: nameless arg {arg}")]
    NamelessArg { prog: String, arg: String },

    #[error("required arguments missing")]
    MissingRequiredParam,

    #[error("getparam: parameter {0} unknown")]
    ParamNotAvailable(String),

    #[error("getparam: called before initparam")]
    NotInitialized,

    #[error("getbparam: {name}={value} not bool")]
    BadBoolParam { name: String, value: String },

    #[error("cannot open file \"{0}\"")]
    FileOpen(String),

    #[error("cannot create file \"{0}\"")]
    FileCreate(String),

    #[error("outputdata: cannot open output file")]
    OutputFileOpen,

    #[error("write failed")]
    WriteFailed,

    #[error("restorestate: fread failed")]
    ReadFailed,

    #[error("in_int: input conversion error")]
    InputIntConversion,

    #[error("in_real: input conversion error")]
    InputFloatConversion,

    #[error("out of memory: {0} bytes")]
    OutOfMemory(usize),

    #[error("cputime: times() call failed")]
    CpuTimeFailed,

    #[error("help")]
    Help,
}

pub fn eprintf(fmt: &str) {
    eprint!("{}", fmt);
    let _ = std::io::stderr().flush();
}
