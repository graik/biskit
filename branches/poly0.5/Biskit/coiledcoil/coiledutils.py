from string import letters,digits

def scoreAlignFun(a,b):
    return abs(ord(a)-ord(b) / max(ord(a),ord(b)))
    
def toStringChain(a,b):
    printable = letters + digits
    
    ar = "";br = ""
    maxp = max(max(a),max(b))
    na = array(a);nb = array(b)
    na = na / maxp;nb = nb / maxp
    maxi = len(printable)-1
    for i in na:
        ar += printable[int(maxi*i)]
    for i in nb:
        br += printable[int(maxi*i)]
    return ar,br
    