import os, sys

try:
    import bldutil
    glob_build = True # scons command launched in RSFSRC
    srcroot = '../..' # cwd is RSFSRC/build/user/rlwalker
    Import('env bindir libdir pkgdir')
    env = env.Clone()
except:
    glob_build = False # scons command launched in the local directory
    srcroot = os.environ.get('RSFSRC', '../..')
    sys.path.append(os.path.join(srcroot,'framework'))
    import bldutil
    env = bldutil.Debug() # Debugging flags for compilers
    bindir = libdir = pkgdir = None
    SConscript(os.path.join(srcroot,'su/lib/SConstruct'))

targets = bldutil.UserSconsTargets()

# C mains
targets.c = '''
pwefd2D
'''

# vse
# awefd2d awefd3d

# Python targets
#targets.py = '''
#'''

dynlib = env.get('DYNLIB','')

env.Prepend(LIBS=[dynlib+'su'])

targets.build_all(env, glob_build, srcroot, bindir, libdir, pkgdir)
