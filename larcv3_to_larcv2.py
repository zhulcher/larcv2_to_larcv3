import os
import ctypes


dir_path = os.path.dirname(os.path.realpath(__file__))
lib = ctypes.cdll.LoadLibrary(dir_path  + "/" + "larcv3_to_larcv2.so")
class larcv3_to_larcv2(object):
    def __init__(self):
        lib.larcv3_to_larcv2_new.argtypes = []
        lib.larcv3_to_larcv2_new.restype = ctypes.c_void_p
        lib.larcv3_to_larcv2_add_in_file.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        lib.larcv3_to_larcv2_add_in_file.restype = ctypes.c_void_p
        lib.larcv3_to_larcv2_set_out_file.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        lib.larcv3_to_larcv2_set_out_file.restype = ctypes.c_void_p
        lib.larcv3_to_larcv2_initialize.argtypes = [ctypes.c_void_p]
        lib.larcv3_to_larcv2_initialize.restype = ctypes.c_void_p
        lib.larcv3_to_larcv2_convert.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
        lib.larcv3_to_larcv2_convert.restype = ctypes.c_void_p
        self.obj = lib.larcv3_to_larcv2_new()

    def add_in_file(self, in_file):
        # in_file = in_file.encode('utf-8')
        # # in_file = bytes(in_file, encoding="ascii")
        # print(in_file)
        # print(type(in_file))
        lib.larcv3_to_larcv2_add_in_file(self.obj, in_file.encode('utf-8'))
    
    def set_out_file(self, out_file):
        lib.larcv3_to_larcv2_set_out_file(self.obj, out_file.encode('utf-8'))

    def initialize(self):
        lib.larcv3_to_larcv2_initialize(self.obj)

    def convert(self, n_events, n_skip):
        lib.larcv3_to_larcv2_convert(self.obj, n_events, n_skip)


import argparse

def main():
    parser = argparse.ArgumentParser(description='LArCV2 to larcv3 Conversion script')

    parser.add_argument('-il','--input-larcv',required=True,
                        dest='larcv_fin',nargs='+',
                        help='string or list, Input larcv file name[s] (Required)')

    parser.add_argument('-nevents','--num-events',
                        type=int, dest='nevents', default=-1,
                        help='integer, Number of events to process')

    parser.add_argument('-nskip','--num-skip',
                        type=int, dest='nskip', default=0,
                        help='integer, Number of events to skip before processing')

    parser.add_argument('-ol','--output-larcv',default='',
                        type=str, dest='larcv_fout',
                        help='string,  Output larcv file name (optional)')

    args = parser.parse_args()
    print(args)
    converter = larcv3_to_larcv2()

    for file_name in args.larcv_fin:
        converter.add_in_file(file_name)
    
    if args.larcv_fout == "":
        print("No output file specified, using basename and same path as first input.")
        args.larcv_fout = args.larcv_fin[0].replace(".h5", ".root")


    converter.set_out_file(args.larcv_fout)

    converter.initialize()
    converter.convert(args.nevents, args.nskip)


if __name__ == "__main__":
    main()


# # help(conversion_libs)
# # print(conversion_libs)

# # conversion_libs.add_in_file