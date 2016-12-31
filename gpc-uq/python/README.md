# Python gPC UQ package

This is a Python implementation of the Matlab package that shares the same repo (see ```../matlab/README.md```).

It is strongly recommended that you test this package in a [virtual environment](https://virtualenv.pypa.io/en/stable/); e.g. by running the following before you start working with this code (assumes Bash shell)

```
$ pip install --upgrade virtualenv
$ virtualenv venv
$ source venv/bin/activate
```

Once you've set up your environment, you can install the gpc package and its dependencies with 

```
(venv) $ pip install --editable .
```

Try running the tests/examples with:

```
(venv) $ python -m unittest -v
```

Sometimes (if you are like the person writing these instructions) you will need to completely start over with a clean environment.
You can do this by uninstalling the gpc module

```
(venv) $ pip uninstall .
```

or deactivating and deleting your virtual environment 

```
(venv) $ deactivate
$ rm -rf venv
```