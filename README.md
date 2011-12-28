
Dynamic Systems Simulator
=========================


Requirements
============

 - Python 2.6+
 - Scipy 0.99
 - Mathplotlib
 - Ph.D. Physics (optional)

Install
=======

For devs

    $ sudo pip install virtualenv
    $ virtualenv --no-site-packages env
    $ source env/bin/activate
    (env)$ pip install {{dependencias}}
    (env)$ git clone git@github.com:horatz/superharm.git superharm-project
    (env)$ cd superharm-project
    (env) ~/superharm-project $ python setup.py install
    (env) ~/superharm-project $ Superharm
    (env) ~/superharm-project $ ./Superharm
    ... varios cambios despues

    (env) ~/superharm-project $ git add superharm
    (env) ~/superharm-project $ git commit -m "Added some strange and awesome feature"
    (env) ~/superharm-project $ git push 

Distributing
============

For windows

    $ python setup.py bdist --formats wininst

