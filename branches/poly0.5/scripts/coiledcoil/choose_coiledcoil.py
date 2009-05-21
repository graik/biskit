from Biskit.polysys import coiledcoil


if __name__ == '__main__':
    
    # Check if arguments are ok
    if len(sys.argv)<4:
        print "\nUsage: exe arg1 \n"
        print "arg1->  File containing data for comparison."
        exit( -1 )
    
    #Open data files...
    try:
        ccdata = open(sys.argv[1],"r")
    except (IOError,EOFError):
        print "-Error opening "+sys.argv[1]+": file doesn't exist."
        exit( -1 )