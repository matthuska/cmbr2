# cmbr

Miscellaneous scripts for working with Bioinformatics data in R.
They're mostly for working with ChIP-seq data but feel free to add R
functions that you think other people might want to use.

## Coding style

Try to follow the coding style of the other code in the repository. If
you're not sure what to do, I'd recommend using Hadley Wickham's
coding style:

https://github.com/hadley/devtools/wiki/Style

At the very least try to do two things:

* indent code with two spaces
* use <- as an assignment operator
* document the function so that other people have some idea what it is
  supposed to do, the arguments that it takes and what the function
  returns

## Adding cmbr as a submodule of an existing git project:

```bash
cd <existing git project root>
git submodule add git@github.com:matthuska/cmbr.git cmbr
git submodule update --init
git commit ./cmbr -m "Added submodule as ./cmbr"
```

## Updating your repository with other people's changes

```bash
cd ./cmbr
git pull --rebase
cd ..
git commit ./cmbr -m "Updated submodule cmbr"
```

## Making changes and pushing them back

```bash
cd ./cmbr
<edit files>
git commit -a -m "<describe your changes>"
git pull --rebase
git push
cd ..
git commit ./cmbr -m "Updated submodule cmbr"
```
