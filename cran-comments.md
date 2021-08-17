## Resubmission
This is a resubmission. In this version:
* An issue was addressed which caused an ERROR when the package was built on 
some platforms (r-patched-solaris-x86).

## Test environments
* R-hub windows-x86_64-devel (r-devel)
* R-hub ubuntu-gcc-devel (r-devel)
* R-hub debian-gcc-devel (r-devel)
* R-hub macos-highsierra-release-cran (r-release)
* R-hub solaris-x86-patched-ods (r-release)

## R CMD check results
There were no ERRORS or WARNINGS.

There was 1 NOTE for all platforms:

* Found the following (possibly) invalid DOIs:
  DOI: 10.2307/2533197
    From: DESCRIPTION
    Status: Forbidden
    Message: 403

This DOI is correct and works.
