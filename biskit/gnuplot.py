## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

# Simple Gnuplot interface.
#
# Written by Konrad Hinsen <hinsen@ibs.ibs.fr>
# last revision: 1998-1-19
#
# Caution: If you use a Gnuplot version earlier than 3.6beta,
# every call for a screen display creates another gnuplot
# process; these processes are never closed. There seems to be no
# other way to make gnuplot behave as it should.

# Modified by Niklas Blomberg, EMBL to include box and scatterplots

# Modified by Raik Gruenberg, IP:
#   allow importing even if gnuplot is not installed,
#   check gnuplot.installed == 1, to be sure the program runs

"""
Simple Gnuplot interface.
"""

import os, string, tempfile

installed = 0
old_version = 0

#
# Test if Gnuplot is new enough to know the option -persist
#
try:
    filename = tempfile.mktemp()
    file = open(filename, 'w')
    file.write('\n')
    file.close()
    gnuplot = os.popen('gnuplot -persist ' + filename + ' 2>&1', 'r')
    response = gnuplot.readlines()
    gnuplot.close()
    os.unlink(filename)
    old_version = response and string.index(response[0], '-persist') >= 0
    installed = 1
except:
    pass

#
# List in which pipes to Gnuplot processes are stored to prevent them from
# being closed at the end of the plot function.
#
_gnuplot_pipes = []

#
# Check if data object is a sequence
#
def _isSequence(object):
    n = -1
    try: n = len(object)
    except: pass
    return n >= 0

#
# Generate the Gnuplot data files and command string for standard plots
#
def _plotData(data):
    plotlist = []
    filelist = []
    for set in data:
        filename = tempfile.mktemp()
        file = open(filename, 'w')
        is_sequence = _isSequence(set[0])
        for point in set:
            if is_sequence:
                for coordinate in point:
                    file.write(repr(coordinate) + ' ')
            else:
                file.write(repr(point))
            file.write('\n')
        file.close()
        if is_sequence:
            plotlist.append((filename, len(set[0])))
        else:
            plotlist.append((filename, 1))
        filelist.append(filename)
    command = 'plot '
    for item in plotlist:
        filename, n = item
        if n == 1:
            command = command + '"' + filename + '"  notitle w l, '
        else:
            for i in range(n-1):
                command = command + '"' + filename + \
                        '"  using 1:' + repr(i+2) + ' notitle w l, '
    command = command[:-2] + '\n'
    return command, filelist

#
# plot of multiple lists with labels
#
def _plotWithLabels(data):
    plotlist = []
    filelist = []
    if not type(data).__name__ == 'dictionary':
        data = zip(data, len(data)*[''])
    else:
        data = list(data.items())
    for set, key in data:
        filename = tempfile.mktemp()
        file = open(filename, 'w')
        is_sequence = _isSequence(set[0])
        for point in set:
            if is_sequence:
                for coordinate in point:
                    file.write(repr(coordinate) + ' ')
            else:
                file.write(repr(point))
            file.write('\n')
        file.close()
        if is_sequence:
            plotlist.append((key, filename, len(set[0])))
        else:
            plotlist.append((key, filename, 1))
        filelist.append(filename)
    command = 'plot '
    for item in plotlist:
        key, filename, n = item
        if n == 1:
            command = command + '"' + filename + '"  title "%s" w l, ' %key
        else:
            for i in range(n-1):
                command = command + '"' + filename + \
                        '"  using 1:' + repr(i+2) + ' title "%s" w l, ' %key
    command = command[:-2] + '\n'
    return command, filelist

###########################################################################
# Barplot of data
###########################################################################

def _barGraphPlotData(data):
    """
    _barGraphPlotData(data):
    """
    plotlist = []
    filelist = []
    for set in data:
        filename = tempfile.mktemp()
        file = open(filename, 'w')
        is_sequence = _isSequence(set[0])
        for point in set:
            if is_sequence:
                for coordinate in point:
                    file.write(repr(coordinate) + ' ')
            else:
                file.write(repr(point))
            file.write('\n')
        file.close()
        if is_sequence:
            plotlist.append((filename, len(set[0])))
        else:
            plotlist.append((filename, 1))
        filelist.append(filename)
    command = 'plot '
    for item in plotlist:
        filename, n = item
        if n == 1:
            command = command + '"' + filename + '" notitle with boxes , '
        else:
            for i in range(n-1):
                command = command + '"' + filename + \
                        '"  using 1:' + repr(i+2) + 'notitle with boxes , '
    command = command[:-2] + '\n'
    #print command
    return command, filelist

###########################################################################
# Scatterplots
###########################################################################

def _scatterData(data, marker = 'points'):
    plotlist = []
    filelist = []
    for set in data:
        filename = tempfile.mktemp()
        file = open(filename, 'w')
        is_sequence = _isSequence(set[0])
        for point in set:
            if is_sequence:
                for coordinate in point:
                    file.write(repr(coordinate) + ' ')
            else:
                file.write(repr(point))
            file.write('\n')
        file.close()
        if is_sequence:
            plotlist.append((filename, len(set[0])))
        else:
            plotlist.append((filename, 1))
        filelist.append(filename)
    command = 'plot '
    for item in plotlist:
        filename, n = item
        if n == 1:
            command = command + '"' + filename + '"  notitle  with %s, ' \
                    %marker 
        else:
            for i in range(n-1):
                command = command + '"' + filename + \
                        '"  using 1:' + repr(i+2) + ' notitle  with %s, ' \
                        %marker
    command = command[:-2] + '\n'
    return command, filelist


# additional

def _scatterData3D(data):
    plotlist = []
    filelist = []
    for set in data:
        filename = tempfile.mktemp()
        file = open(filename, 'w')
        is_sequence = _isSequence(set[0])
        for point in set:
            if is_sequence:
                for coordinate in point:
                    file.write(repr(coordinate) + ' ')
            else:
                file.write(repr(point))
            file.write('\n')
        file.close()
        if is_sequence:
            plotlist.append((filename, len(set[0])))
        else:
            plotlist.append((filename, 1))
        filelist.append(filename)
    command = 'splot '
    for item in plotlist:
        filename, n = item
        if n == 1:
            command = command + '"' + filename + '"  notitle  with points, '
        else:
            for i in range(n-1):
                command = command + '"' + filename + \
                        '"  using 1:' + repr(i+2) + ' notitle  with points, '
    command = command[:-2] + '\n'
    return command, filelist



#
# Generate the Gnuplot data files and command string for parallel-axes plots
#
def _parallelAxesPlotData(data, origin):
    naxes = len(data[0])
    filename = tempfile.mktemp()
    filelist = [filename]
    file = open(filename, 'w')
    lower = data[0][0]
    upper = lower
    for point in data:
        for i in range(naxes):
            value = point[i]
            file.write('%d %g\n' % (i+origin, value))
            lower = min(lower, value)
            upper = max(upper, value)
        file.write('\n')
    margin = 0.05*(upper-lower)
    for i in range(0, naxes):
        file.write('\n')
        file.write('%d %g\n' % (i+origin, lower-margin))
        file.write('%d %g\n' % (i+origin, upper+margin))
    file.close()
    command = 'plot ' + '[%d:%d] [%g:%g]' % (origin, origin+naxes-1,
                                             lower-margin, upper+margin) + \
            '"' + filename + '" index 0 notitle w l, "' + \
            filename + '" index 1 notitle w l lt -1\n'
    return command, filelist

#
# Execute a Gnuplot command
#
def _execute(command, filelist, keywords):
    if 'file' in keywords:
        filename = tempfile.mktemp()
        file = open(filename, 'w')
        file.write('set terminal postscript\n')
        file.write('set output "' + keywords['file'] + '"\n')
        file.write(command)
        file.close()
        filelist.append(filename)
        os.system('gnuplot ' + filename)
    else:
        if old_version:
            gnuplot = os.popen('gnuplot 1> /dev/null 2>&1', 'w')
            gnuplot.write('set terminal x11\n')
            gnuplot.write(command)
            gnuplot.flush()
            _gnuplot_pipes.append(gnuplot)
            os.system('sleep 2s')
        else:
            gnuplot = os.popen('gnuplot -persist 1> /dev/null 2>&1', 'w')
            gnuplot.write('set terminal x11\n')
            gnuplot.write(command)
            gnuplot.write('quit\n')
            gnuplot.close()
    for file in filelist:
        os.unlink(file)



#
# Generate a plot
#
def plot(*data, **keywords):
    command, filelist = _plotData(data)
    _execute(command, filelist, keywords)


def plotWithLabels(data, **keywords):
    command, filelist = _plotWithLabels(data)
    _execute(command, filelist, keywords)
#
#Generate Boxplot
#
def boxPlot(*data, **keywords):
    """
    boxPlot
    """
    command, filelist =  _barGraphPlotData(data)
    _execute(command, filelist, keywords)
    return
#
#Scatterplot
#
def scatter(*data, **keywords):
    try:
        marker = keywords['marker']
    except:
        marker = 'points'
    command, filelist = _scatterData(data, marker = marker)
    _execute(command, filelist, keywords)


#
# Generate a parallel-axes plot
#
def parallelAxesPlot(data, **keywords):
    try:
        origin = keywords['origin']
    except KeyError:
        origin = 0
    command, filelist = _parallelAxesPlotData(data, origin)
    _execute(command, filelist, keywords)


#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test case"""

    def prepare(self):
        self.fout = tempfile.mktemp('ps','testgnuplot_')

    def cleanUp(self):
        import biskit.tools as T
        T.tryRemove( self.fout )

    def test_plot2ps(self):
        """gnuplot.plot to file test"""
        plot([1, 5, 3, 4], file = self.fout)
        if self.local:
            print('plot written to ', self.fout)

    def test_scatter(self):
        """gnuplot.scatter test (interactive only)"""
        from numpy.random.mtrand import poisson
        if self.local:
            self.p = scatter( poisson(50,(1000,2))  )

    def test_plot( self ):
        """gnuplot.plot test"""
        # List of (x, y) pairs
        # plot([(0.,1),(1.,5),(2.,3),(3.,4)])
        # plot( zip( range(10), range(10) ) )

        # Two plots; each given by a 2d array
        import biskit.core.oldnumeric as N0
        x = N0.arange(10)
        y1 = x**2
        y2 = (10-x)**2
        plot( N0.transpose(N0.array([x, y1])), N0.transpose(N0.array([x, y2])))

    def test_parallelAxesPlot(self):
        """gnuplot.parallelAxesPlot test (interactive only)"""
        if self.local:
            data = [[0., 1., 0.], [1., -1., 0.], [0.5, 0., 1.]]
            parallelAxesPlot(data)


if __name__ == '__main__':

    BT.localTest()
