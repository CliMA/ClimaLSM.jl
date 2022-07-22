var documenterSearchIndex = {"docs":
[{"location":"#ClimaLSM.jl","page":"Home","title":"ClimaLSM.jl","text":"","category":"section"},{"location":"Contributing/#Contributing","page":"Contribution guide","title":"Contributing","text":"","category":"section"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"Thank you for contributing to ClimaLSM! We encourage Pull Requests (PRs). Please do not hesitate to ask questions.","category":"page"},{"location":"Contributing/#Some-useful-tips","page":"Contribution guide","title":"Some useful tips","text":"","category":"section"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"When you start working on a new feature branch, make sure you start from master by running: git checkout master.\nMake sure you add tests for your code in test/ and appropriate documentation in the code and/or in docs/. All exported functions and structs must be documented.\nWhen your PR is ready for review, clean up your commit history by squashing and make sure your code is current with ClimateMachine master by rebasing.","category":"page"},{"location":"Contributing/#Continuous-integration","page":"Contribution guide","title":"Continuous integration","text":"","category":"section"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"After rebasing your branch, you can ask for review. Fill out the template and provide a clear summary of what your PR does. When a PR is created or updated, a set of automated tests are run on the PR in our continuous integration (CI) system.","category":"page"},{"location":"Contributing/#Automated-testing","page":"Contribution guide","title":"Automated testing","text":"","category":"section"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"Currently a number of checks are run per commit for a given PR.","category":"page"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"JuliaFormatter checks if the PR is formatted with .dev/climaformat.jl.\nDocumentation rebuilds the documentation for the PR and checks if the docs are consistent and generate valid output.\nTests runs the file test/runtests.jl,  using Pkg.test(). These are a mix of unit tests and fast integration tests.","category":"page"},{"location":"Contributing/","page":"Contribution guide","title":"Contribution guide","text":"We use bors to manage merging PR's in the the ClimaLSM repo. If you're a collaborator and have the necessary permissions, you can type bors try in a comment on a PR to have integration test suite run on that PR, or bors r+ to try and merge the code.  Bors ensures that all integration tests for a given PR always pass before merging into master.","category":"page"}]
}