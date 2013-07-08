\name{boom}
\alias{boom}
\title{Segfault Test Function boom}
\usage{
  boom()
}
\value{
  R will crash with a segfault but you will see a
  traceback. ONLY run in console.
}
\description{
  test function produces a C++ segfault. Calls function baz
  which allocates a wrong pointer. If you compiled this
  code with \code{CXXFLAGS=-g3 -rdynamic} the installed
  function \code{handler} will print a traceback of the
  stack that contains the name of the offending function.
  Without this compiler flag, you miss the function name.
  you should place a call to \code{signal(SIGSEGV,
  handler);} at the beginning of each function you want to
  check for segfaults.
}

