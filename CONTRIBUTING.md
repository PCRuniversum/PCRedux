# How to contribute

You would like to contribute? We are happy to see you here.

## Testing

`PCRedux` is tested with the [testthat](https://cran.r-project.org/package=testthat) package. Please write new tests (where applicable) for any code you create.

Our package contains a lot of datasets that were used to test the package. However, qPCR amplification curves are highly diverse. This testing with your data is central to improve the quality of the software.

## Submitting Changes

Once you're ready to submit code or data, please create a pull request that clearly outlines the changes you've made. Additionally, when committing code please use
clear commit messages describing the feature being committed. This helps keep track of changes to the project and makes reviews much easier!

The package used a MIT license.

## Coding Conventions

- We use mainly underscore_separated code style
- The package is written in S3
- If you think the code isn't simple enough to understand alone, write comments!
- Please use the [styler](https://cran.r-project.org/package=styler) to improve the readbility of the code
