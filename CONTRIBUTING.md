
# Contribution guidelines

First of all, thank you very much for looking into how to contribute!

We are delighted to have everyone involved in the project and would like you to feel welcome!

## What is a contribution

I am defining contributing in a very broad way, so any form of input that helps the project
in the future is a contribution.

There are many ways to help the project:
1. Reporting problems
2. Suggesting improvement
3. Discussing design decisions
4. Adding training data
5. Testing the software
6. Adding code
7. ...

## Contribution workflow

In general it is a great idea to start the contribution with an issue.
In the issue the general idea is to specify what the expected outcome would be (start a discussion, fix an error, offer help writting a tutorial ...). After this we will reply and get on the same page on the best way to get it done!

No issue is too small, even an issue saying "I do not understand how to use this project" is a valuable source of feedback.


## Code guidelines

### Formatting

We have several tools that check for the code standards, most of them are automated using pre-commit.
If you want to run all the tools on all the files use this:

```shell
pipx install pre-commit
pre-commit install
pre-commit --all-files
```

### Testing

We intend to maintain a fairly high code coverage (Ideally 100% but there are exceptions).
The coverage can be checked using:

```shell
pytest --cov-report html
```

This will generate a directory named `htmlcov` which contains the coverage report.

Please note that 100% coverage is not a perfect proxy for a piece of software without bugs.
In general, adding functionality needs to be in conjunction with tests that check
for that functionality explicitly.

### Versioning

We use semantic versioning

```shell
bumpver update --minor
```

### Publishing

Publishing is handled using github actions but the underlying
commands for it are these (you should never have to run this manually)

```shell
rm -rf dist/*
python -m build
twine upload dist/*
```

### Maintaining dependencies

asdad
