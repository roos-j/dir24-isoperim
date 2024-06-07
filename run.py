from argparse import ArgumentParser
from dir24isoperim import verify_all, b0, b1, c0, init_prec, Output
from flint import arb, ctx

if __name__ == "__main__":
    parser = ArgumentParser(description="Verify estimates in DIR24.")
    parser.add_argument("--beta", type=str, default=b0, dest="beta",
                        help="Value for beta0; should be in [0.5, %.3f] (default: %f)"%(float(b1),float(b0)))
    parser.add_argument("--c", type=str, default=c0, dest="c",
                        help="Value for c0; should be in [0,1] (default: %f)"%float(c0))
    parser.add_argument("--prec", type=int, default=ctx.prec, dest="prec",
                        help="Working precision in bits (default: %d)"%ctx.prec)
    parser.add_argument("--filename", type=str, default="", dest="filename", 
                        help="Write partition data to a file.")
    args = parser.parse_args()
    
    #ctx.prec = args.prec
    init_prec(args.prec)
    print("Working precision: %d"%ctx.prec)
    beta = arb(args.beta)
    c = arb(args.c)
    print("beta0 = %s"%beta)
    print("c0 = %s"%c)
    if args.filename: 
        if Output.get_instance().open(args.filename):
            print("Partition data will be written to '%s'"%args.filename)

    verify_all(beta, c)

    Output.get_instance().close() 
