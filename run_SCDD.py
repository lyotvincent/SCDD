import argparse

def shellAccept():
    try:
        parser = argparse.ArgumentParser(description='SCDD--Single Cell Diffusion and Denoising.')
        parser.add_argument("-n", "--name", type=str, default="counts",
                            help="the name of experiments or dataset to identify the imputed results")
        parser.add_argument("-r", "--raw", type=str, default=None, help="the raw data path, except format: `.tsv` and \
                        `.h5ad`, make sure that your delimiter in .tsv file is `\t` and your raw data is .X in .h5ad file.")
        parser.add_argument("-tr", "--trans", type=bool, default=True, help="if raw data layout is: columns are cells \
                        and rows are genes, set this param `True`, else `False`.")
        parser.add_argument("-i", "--id", type=int, help="the id number to identify different results from the same methods, \
                            default 0.", default=0)
        parser.add_argument("-f", "--format", type=str, help="the format of input raw data file, support `tsv` and `h5ad`, \
                            default None, this will automatically identify the format by the suffix of the raw file.", default=None)
        parser.add_argument("-me", "--max-epoch", type=int, help="the max epoch of Denoising, default 500", default=500)
        parser.add_argument("-nb", "--neighbors", type=int, help="the neighbors of each cell, default 20", default=20)
        parser.add_argument("-t", "--threshold", type=float, help="the threshold of whether the point is a drop-out point, the \
                            drop-out rate higher than this value will be treated as a drop-out point, default 0.2.", default=0.2)
        parser.add_argument("-nm", "--neighbor-method", type=str, help="the method to generate concensus matrix and neighbors \
                            relations, default `SC3`. if the number of cells is more than 10000, it will automatically turn to `SNN`.",
                            default="SC3")
        parser.add_argument("-b", "--batch-size", type=int, help="the batch_size of each input in the nerual network when \
                            denoising, default 5000, which means if cell number are less than 5000, the total batch will be 1.",
                            default=5000)
        return parser
    except Exception as e:
        print(e)
        exit(1)


if __name__ == '__main__':
    parser = shellAccept()
    args = vars(parser.parse_args())
    from SCDD_api import *
    e = SCDD(name=args['name'],
             raw=args['raw'],
             id=args['id'],
             Tran=args['trans'],
             format=args['format'],
             batch_size=args['batch_size'],
             structure="AnnData",
             neighbor_method=args['neighbor_method'],
             max_epoch=args['max_epoch'],
             neighbors=args['neighbors'],
             threshold=args['threshold'],)
    e.run(store=False)






