Developer Info
==============

To get started developing you obviously need git installed. You can make
changes in your local directory, and if you make a change that you think
is substantial and general enough that others would be interested you
can "push" the change using git (see the `git
tutorial <http://git-scm.com/docs/gittutorial>`__). If you are going to
submit a change, please update the version number in
catmap/\ **init**.py (see **version** == x.x.x). The first number in the
version is used for major changes, and a 0 represents code which is
still not completely stable. The second number represents moderate
changes (significant new functionality added), and the final number
corresponds to the git submission number. The final number is the most
annoying to update, since you need to change it every time before you
submit. In order to make this slightly easier you can use the following
script:

.. code:: python

    from subprocess import Popen, PIPE

    output = Popen(['git', 'rev-list', 'HEAD','--count'],stdout=PIPE)
    rev = output.stdout.read().strip()
    init = open('catmap/__init__.py')
    init_txt = init.read()
    init.close()
    new_txt = []
    for li in init_txt.split('\n'):
        if '__version__' in li:
            oldversion = li.rsplit('.',1)[0]
            newversion = oldversion + '.' + rev + '"'
            new_txt.append(newversion)
        else:
            new_txt.append(li)

    new_init = '\n'.join(new_txt)
    init = open('catmap/__init__.py','w')
    init.write(new_init)
    init.close()

    message = raw_input('Submit message:')
    os.system('git commit -am "'+message+'"')
    os.system('git push')

Place the script in the base directory of the project and run it with
python. If anyone knows a better way of keeping the version number
synchronized I am open to suggestions.
