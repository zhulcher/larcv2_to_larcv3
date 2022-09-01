import filecmp
import subprocess

  
f1 = "test2_3.root"
f2 = "test2_3.h5"

f3 = "test3_2.h5"
f4 = "test3_2.root"

subprocess.run(["python3", "larcv2_to_larcv3.py","-il",f1])
subprocess.run(["cp",f2,f3])
subprocess.run(["python3", "larcv3_to_larcv2.py","-il",f3])

def test_3to2():
    assert filecmp.cmp(f1, f4, shallow=False)


subprocess.run(["python3", "larcv3_to_larcv2.py","-il",f2])
subprocess.run(["cp",f1,f4])
subprocess.run(["python3", "larcv2_to_larcv3.py","-il",f4])

def test_3to2():
    assert filecmp.cmp(f2, f3, shallow=False)