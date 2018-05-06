import datetime
import itertools
import random

class clsQCreator(object):

    """
    Holds discharge and date time series.  Differs from clsTimeSeries in that it
    runs in pypy
    """

    def __init__(self, DateList, QList):

        """
        Attributes:

        -Datelist--[tuple] (List of date tuples; will be converted to datetime objects)
        -Qlist [float] (List of discharge values)
        -numH2OYears--int (Number of complete water years in record)
        """

        self.Dates = map(lambda x: datetime.datetime.strptime(x, '%Y,%m,%d').date(), DateList)
        self.QList = map(lambda x: float(x), QList)
        self.StartH2OYear = ""
        self.numH2Oyears = ""
        self.FindNumH2OYears()

    def FindNumH2OYears(self):
        """
        Outputs the number of complete water years in the object attribute's record
        """
        # Find number of years in record
        StartH2OYear = self.Dates[0].year + 1
        EndH2OYear = self.Dates[-1].year

        if datetime.date(EndH2OYear, 9, 30) not in self.Dates: # Only want to randomize complete water years
            EndH2OYear = EndH2OYear - 1
        if datetime.date(StartH2OYear - 1, 10, 1) not in self.Dates:
            StartH2OYear = StartH2OYear + 1

        self.StartH2OYear =  StartH2OYear   
        self.numH2Oyears = EndH2OYear - StartH2OYear + 1

    def RandomizeWaterYears(self, MaxPermutations = 'All'):

        """
        Creates discharge records in which the order of daily discharges within
        individual water years is preserved, but the order of the years is randomized

        Right now, this function only supports creating the max number of permutations.

        Attributes:

        -StartH20Year--int (Water year in which randomization starts--if 'Beginning,
            all complete years will be randomized)
        -MaxPermutations--int (Number of random records to create.  If 'All', all
            permutations will be created.)
        """
        # Find last water year
        LastYear = self.Dates[-1].year
        EndH2OYear = LastYear

        if datetime.date(LastYear, 9, 30) not in self.Dates: # Only want to randomize complete water years
            EndH2OYear = LastYear - 1

        numyears = EndH2OYear - self.StartH2OYear

        # Create all the permutations of years and create Q records from them 
        permutations = list(itertools.permutations(range(numyears)))

        return permutations

    def RandomizeRecord(self, order, StartH2OYear = 'Beginning'):

        """
        Takes a permutation and changes the order of discharges in the QList
        attribute to reflect the order in the permutation.

        Attributes:

        order--(int) (Tuple of the random order of years)
        StartH2OYear--int (Year to start randomizing)
        """

        LastYear = self.Dates[-1].year
        EndH2OYear = LastYear
        QYearDict = {}

        if datetime.date(LastYear, 9, 30) not in self.Dates: # Only want to randomize complete water years
            EndH2OYear = LastYear - 1

        numyears = EndH2OYear - StartH2OYear + 1
        startindex = 0
        endindex = self.Dates.index(datetime.date(EndH2OYear, 10, 1))

        # Give each discharge a key
        for i in range(numyears):
            index1 = self.Dates.index(datetime.date(StartH2OYear + i - 1, 10, 1))
            index2 = self.Dates.index(datetime.date(StartH2OYear + i, 10, 1))

            if i == 0:
                startindex = index1

            QYearDict[i] = self.QList[index1:index2]

        # Create the new discharge record
        NewQList = self.QList[0:startindex]
        for num in list(order):
            NewQList = NewQList + QYearDict[num]
        NewQList = NewQList + self.QList[endindex:]

        self.QList = NewQList

    def FabricateRandomOrder(self, numyears, numpermutations):
        """
        Creates random orders of discharges with of a user-defined number of years.  Draws from all
        complete water years in the object 'Qlist' field.

        Attributes:
        -numyears (int)--number of years that the output record will be.
        -numpermutations (int)--number of permutations to export
        """
        permutationdict = {}
        for i in range(numpermutations):
            randomlist = []
            for j in range(numyears):
                randomlist.append(random.randint(0, self.numH2Oyears-1))
            permutationdict[i]=randomlist
        
        return permutationdict

    def FabricateRandomRecord(self, permutation):
        """
        Creates a random record using the annual hydrographs in the 'Qlist' field, given a user-
        provided list of indexes.

        Attributes:
        -permutation ([int])--list of random integers, which cannot be larger than number of years
        in record
        """
        RandomQList = []
        
        # Break up QList by year
        QYearDict = {}

        # Give each discharge a key
        for i in range(self.numH2Oyears):
            index1 = self.Dates.index(datetime.date(self.StartH2OYear + i - 1, 10, 1))
            index2 = self.Dates.index(datetime.date(self.StartH2OYear + i, 10, 1))

            QYearDict[i] = self.QList[index1:index2]

        # Create discharge record using permutations
        for i in permutation:
            RandomQList = RandomQList + QYearDict[i]

        return RandomQList

        












        
        
