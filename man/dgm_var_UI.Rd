% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_modules.R
\name{dgm_var_UI}
\alias{dgm_var_UI}
\title{dgm_var_UI}
\usage{
dgm_var_UI(id)
}
\arguments{
\item{id, }{character used to specify namespace, see \code{shiny::\link[shiny]{NS}}}
}
\value{
a \code{shiny::\link[shiny]{tagList}} containing UI elements
}
\description{
Variables UI control selection for persistence diagrams.
}
\details{
\code{shiny::\link[shiny]{tagList}} returned includes a selectizeInput, a numericInput, and
a dynamically generated uiOutput that becomes a range slider.
}
