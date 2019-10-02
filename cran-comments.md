## Test environments
* local OS X install, R 3.1.2
* ubuntu 16.04 LTS - R 3.6.1
* local OS X install, R 3.1.2
* Rhub Ubuntu Linux 16.04 LTS

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility

  Version contains large components (0.0.0.9000)

## R submission 

### v0.0.2

Pls use a sensible version number and resubmit. then we will inspect manually.

The Description field is intended to be a (one paragraph) description
of what the package does and why it may be useful.
Please elaborate and explain what entrezID is.

### v0.0.3

If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.

You write information messages to the console that cannot be easily suppressed.
Instead of print()/cat() rather use message()/warning()  or if(verbose)cat(..) if you really have to write text to the console.
(except for print() and summary() functions)

Please fix and resubmit, and document what was changed in the submission comments.


### v0.0.5
The Description field should not start with the package name,
  'This package' or similar.
The Description field contains
  <https://doi.org/10.1093/nar/gkv007> to find the correspondence between
Please write DOIs as <doi:10.prefix/suffix>.

Please fix and resubmit. 

