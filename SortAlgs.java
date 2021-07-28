/**
 * @author      Bogdan Kovch
 * @class       CS 350
 * @assignment  final project
 * @file        SortAlgs.java
 * @date        May 29, 2013 - June 6, 2013
 */

/**
 * This program analyzes one of the sorting algorithms specified by the user. The program
 * automatically generates a series of arrays of integers (sequences) for every test case and
 * feeds them into the algorithm. Every algorithm has built in routines for gathering statistics
 * about algorithm's work while processing a particular sequence of integers. Statistics for
 * different individual sequences is combined into the statistics of an individual test case. Then
 * every test case statistics is passed to a procedure which matches this statistics to various
 * mathematical models of efficiency classes and computes percentage how well this data matches
 * every model. The efficiency model which matched the given statistics with the highest percentage
 * is then selected to describe the behavior of the given test case efficiency for the given
 * algorithm. For more details, please see the project report and comments within code.
 * 
 * Compilation notes: use flags -Xss1024m and -Xmx4096m to deal with possible memory limitations.
 */

package sortalgs;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

public class SortAlgs {

    private static Random rand;           // random number generator for random sequences
    private static int base = 10;         // base for sequence sizes in base^i for 0 < i
    private static int exp_incr = 1;      // exponent increase of i in base^i
    private static final double MAX_TIME = Math.pow(10, 10); // 10^9 ns = 1 s
    private static final int LENGTH_MIN = 1000; // minimum length required for an array to affect
                                                // efficiency class match.

    public enum Order {RANDOM, DESCENDING, ASCENDING} // order of elements in a generated sequence
    public enum Algorithm {INSERTION, QUICK, MERGE}   // sorting algorithms implemented
    public enum EfficiencyClass {NONE, CONSTANT,       // possible efficiency classes
                    LOGARITHMIC, LINEAR, LINEARLOGARITHMIC, QUADRATIC, EXPONENTIAL}
    
    // method:  program entry
    // input:   at least one command line argument
    public static void main(String[] args) {
        rand = new Random((System.currentTimeMillis()));    
        if(args.length < 1) {       // at least one argument is required
            System.out.println("command-line argument expected: 1. 'i', 'm', or 'q'");
            System.out.println("                                2. array length base");
            return;                 // no arguments given, terminate
        }
        if(args.length >= 2)        // array length base for [length = base^exponent]
            base = Integer.decode(args[1]);
        if(args.length >= 3)        // array length exponent increment for [length = base^exponent]
            exp_incr = Integer.decode(args[2]);
                
        Algorithm algorithm;        
        switch(args[0]){            // parse command line argument to define algorithm for testing
            case "i":               // "i" stands for insertion sort
                algorithm = Algorithm.INSERTION; break;
            case "q":               // "q" stands for quick sort
                algorithm = Algorithm.QUICK; break;
            case "m":               // "m" stands for merge sort
                algorithm = Algorithm.MERGE; break;
            default:                // unknown argument provided
                System.out.println("invalid 1st argument - use 'i', 'm', or 'q'.");
                return;
        }
        doAlgorithmAnalysis(algorithm); // do the analysis of the selected algorithm
    }

    ////////////////////////////////////// Analysis Methods //////////////////////////////////////
    
    // method:  manages the analysis process of the algorithm provided
    //          generates test cases one for each sequence order type
    //          generates a series of sequences for every test case
    //          feeds each sequence into a specified sorting algorithm
    //          double checks whether the algorithm actually sorted the sequence
    //          gathers algorithm performance statistics for every sequence in every test case
    //          analises and defines efficiency class for every test case for the algorithm
    //          displays statistics and analysis results
    // input:   sorting algorithm to be analysed
    public static void doAlgorithmAnalysis(Algorithm algorithm) {
        for(Order order: Order.values()) {  // generate test cases for every sequence order type
            System.out.printf("\n%s sort; input order: %s\n", // print desciption of
                    algorithm.toString().toLowerCase(),       // algorithm and 
                    order.toString().toLowerCase());          // sequence order type
            List<Statistics> stats = new ArrayList<>();       // stores test case statistics
            double time = System.nanoTime();                // time when this test case started
            try {                                           // handle possible memory exceptions
                for(int i = exp_incr;                       // generate sequences for sorting
                        System.nanoTime()-time < MAX_TIME;  // until time is gone
                        i += exp_incr) {                    // increase next sequence multiple times
                    int[] sequence = generateSequence(order, (int)Math.pow(base, i));
                    Statistics s;                           // contains statistics for this sequence
                    switch(algorithm) {                     // select requested sorting algorithm 
                        case INSERTION:
                            s = insertionSort(sequence); break;
                        case QUICK:
                            s = quickSortInit(sequence); break;
                        case MERGE:
                            s = mergeSortInit(sequence); break;
                        default:
                            return;
                    }
                    if(!sequenceIsSorted(sequence)) // check if the algorithm sorted the sequence
                        System.out.println("Error! Sequence is not sorted! Algorithm failure!" );
                    stats.add(s); // add sequence sorting statistics to test case statistics
                }
            } catch (OutOfMemoryError e) {
                System.out.println("test case terminated due to \"OutOfMemoryError\" exception");
            }
            printStats(stats); // display current test case statistics
            evaluateEfficiencyClass(stats); //evaluate test case efficiency class
        }
    }
    
    // method:  generates a sequence of specified size with elements in given order
    // input:   order (random, descending or ascending) of elements in the sequence to be generated
    //          size of the sequence to be generated
    // output:  sequence (array of integers) generated
    public static int[] generateSequence(Order order, int size) {
        if(size <= 0)                           // require nonnegative size
            return null;
        int[] sequence = new int[size];         // allocate a new array of needed size
        switch(order){
            case RANDOM:
                for(int i = 0; i < size; i++)
                    sequence[i] = rand.nextInt(Integer.MAX_VALUE);
                return sequence;
            case DESCENDING:
                for(int i = 0; i < size; i++)
                    sequence[i] = size - i;
                return sequence;
            case ASCENDING:
                for(int i = 0; i < size; i++)
                    sequence[i] = i;
                return sequence;
            default:                             // unknown order type
                return null;
        }
    }

    // method:  checks whether the given sequence is sorted in ascending order
    // input:   sequence (array of integers) to be checked
    // output:  true if sequence is sorted or false otherwise
    public static boolean sequenceIsSorted(int[] sequence) {
        if(sequence == null || sequence.length <= 1) // validate sequence passed
            return true;
        int length = sequence.length;
        for(int i = 0, j = 1; j < length; i++, j++)  // compare every pair of consecutive elements
            if( sequence[i] > sequence[j] )
                return false;
        return true;
    }
    
    // method:  displays test case statistics in the form of table
    // input:   list containing statistics about processing a particular test case
    // output:  display output visualizing the statistics
    public static void printStats(List<Statistics> stats) {
        System.out.printf("%10s %15s %15s %15s %15s\n",     // table header
                          "length", "comparisons", "memory ops", "total time", "time/element");
        Iterator<Statistics> iterator = stats.iterator();
        while(iterator.hasNext()) {                         // display every entry in the table
            Statistics s = iterator.next();
            if(s == null || s.length <= 0)                  // validate entry and data
                continue;
            long timePerElement = s.time / s.length;        // calculate field "time/element"
            System.out.printf("%10d %15d %15d %12.3f ms %15s\n",
                              s.length, s.comparisons, s.memoryOps,
                              s.time *  Math.pow(10, -6),
                              timePerElement +  " ns" );
        }
    }
    
    // method:  defines efficiency class for processing the test case described in given statistics
    // input:   list containing statistics about processing a particular test case
    // output:  displays to which efficiency class processing of this test case belongs
    public static void evaluateEfficiencyClass(List<Statistics> stats) {
        EfficiencyClass bestEfficiency = EfficiencyClass.NONE;  // efficiency with the best match
        double bestMatch = 0.0;                                 // best match score
        for(EfficiencyClass efficiency: EfficiencyClass.values()) { // try every efficiency class
            double match = matchEfficiencyClass(stats, efficiency);
            if(match > bestMatch) {                                 // store the best match info
                bestMatch = match;
                bestEfficiency = efficiency;
            }
        }
        System.out.printf("%s: %3.1f%%\n",                          // print the best match info
                          bestEfficiency.toString().toLowerCase(),  bestMatch * 100.0);
    }
    
    // method:  computes the score describing how well the test case statistics matches the
    //          efficiency class offered
    // input:   list containg statistics about a test case
    //          efficiency class to which the statistics needs to be matched and evaluated
    // output:  the score within [0.0, 1.0] describing how the statistics mathes the offered
    //          efficiency class: 0.0 - no match, 1.0 - parfect match
    public static double matchEfficiencyClass(List<Statistics> stats, EfficiencyClass effClass) {
        int size = stats.size();
        double lower = Double.MAX_VALUE;    // lower bound constant of Theta(n^2) 
        double upper = Double.MIN_VALUE;    // upper bound constant of Theta(n^2) 
        
        // Find potential min and max constants for lower and upper bounds of Theta(n^2).
        // If all constants in the statistics within the lower and upper bound constants are very
        // similar in terms of the offered efficiency model then the statistics trend matches this
        // efficiency class well.
        for(int i = 0; i < size; i++) {     // ignore statistics of the first smallest sequences
            Statistics s = stats.get(i);    // retrieve an individual statistics entry
            if(s == null)                   // check if entry is valid
                continue;
            double length = (double)s.length;
            if(length < LENGTH_MIN )        // ignore statistics from small arrays
                continue;            
            double x;                       // potential upper/lower bound constant
            switch(effClass) {              // select math model for offered efficiency class
                case CONSTANT:
                    x = s.comparisons; break;
                case LINEAR:    // solve (x * length = s.comparisons) for x
                    x = s.comparisons / length; break; 
                case LOGARITHMIC:// solve (x * log_ length = s.comparisons) for x
                    x = s.comparisons / (Math.log10(length)/Math.log10(2)); break;
                case LINEARLOGARITHMIC:     // solve (x * length log_2 length = s.comparisons) for x
                    x = s.comparisons / (length * Math.log10(length)/Math.log10(2)); break;
                case QUADRATIC: // solve (x * length^2 = s.comparisons) for x
                    x = s.comparisons / (length * length); break;
                default:
                    return 0.0; // unknown efficiency class
            }
            if(x < lower)
                lower = x;
            if(x > upper)
                upper = x;
        }
        if( upper != 0)
            return lower/upper; 
        else
            return 0.0;
    }
    
    ////////////////////////////////// Insertion Sort Algorithm //////////////////////////////////

    // method:  sorts a given sequence using insertion sort algorithm
    // input:   sequence (array of integers) to be sorted
    // output:  statistics describing the process of sorting the sequence given
    //          sequence sorted in ascending order
    public static Statistics insertionSort(int[] sequence) {
        if(sequence == null)                       // check if sequence is valid
            return null;
        int length = sequence.length;
        Statistics s = new Statistics(length);     // object for storing algorithm's statisitcs
        s.time = System.nanoTime();                // take start time
        for(int i = 1; i < length; i++ ) {         // traverse all sequence
            int value = sequence[i];               // currently selected value to be inserted
            int idx = i;                           // index of the currently selected value
            s.comparisons++;                       // gathering statistics
            while(idx > 0 && value < sequence[idx-1]) { // traverse unsorted part of the sequence
                sequence[idx] = sequence[idx - 1]; // shift an element to the right
                idx--;                             // select position to the left
                s.memoryOps++;                     // gathering statistics
                s.comparisons++;                   // gathering statistics
            }
            sequence[idx] = value;                 // insert the value
            s.memoryOps++;                         // gathering statistics
        }
        s.time = System.nanoTime() - s.time;       // take finish time and compute total time
        return s;
    }
    
    //////////////////////////////////// Quick Sort Algorithm ////////////////////////////////////
    
    // method:  prepares statistics before invoking the quicksort algorithm
    // input:   sequence (array of integers) to be sorted
    // output:  statistics describing the process of sorting the sequence given
    //          sequence sorted in ascending order
    public static Statistics quickSortInit(int[] sequence) {
        if(sequence == null)                        // check if sequence is valid
            return null;
        int length = sequence.length;        
        Statistics s = new Statistics(sequence.length);
        s.time = System.nanoTime();                 // take start time
        quickSort(sequence, 0, length-1, s);        // do the actual quicksort
        s.time = System.nanoTime() - s.time;        // take finish time and compute total tome
        return s;
    }
    
    // method:  sorts a subsequence using quicksort algorithm
    // input:   sequence (array of integers)
    //          left and right indices which define the subsequence of sequence
    //          object for storing algorithm's statistics
    // output:  subsequence sorted in ascending order
    //          updated statistics
    public static void quickSort(int[] sequence, int left, int right, Statistics s) {
        if(left < right) {
            int split = hoarePartition(sequence, left, right, s); // sort this partition
            quickSort(sequence, left, split-1, s);                // sort left partition
            quickSort(sequence, split+1, right, s);               // sort right partition
        }
    }
    
    // method:  partitions and sorts a subarray using Hoare's partitioning algorithm
    // input:   sequence (array of integers)
    //          left and right indices which define the subsequence of sequence
    //          object for storing algorithm's statistics
    // output:  pivot index for further splitting into partitions
    //          updated statistics
    public static int hoarePartition(int[] sequence, int left, int right, Statistics s) {
        int pivot, i, j;
        pivot = sequence[left];
        i = left;
        j = right + 1;
        while(true) {
            s.comparisons += 2;
            while(i < right && sequence[++i] < pivot)
                s.comparisons++;
            while(sequence[--j] > pivot)
                s.comparisons++;
            if(i >= j)
                break;
            s.memoryOps++;
            swap(sequence, i, j);
        }
        s.memoryOps++;
        swap(sequence, left, j);
        return j;
    }

    //////////////////////////////////// Merge Sort Algorithm ////////////////////////////////////
    
    // method:  prepares statistics before invoking the mergesort algorithm
    // input:   sequence (array of integers) to be sorted
    // output:  statistics describing the process of sorting the sequence given
    //          sequence sorted in ascending order
    public static Statistics mergeSortInit(int[] sequence) {
        if(sequence == null)
            return null;
        int length = sequence.length;        
        Statistics s = new Statistics(sequence.length);
        s.time = System.nanoTime();                // take start time
        mergeSort(sequence, s);        // do the actual quicksort
        s.time = System.nanoTime() - s.time;       // take finish time and compute total tome
        return s;
    }
    
    // method:  sorts a sequence using recursive mergesort algorithm
    // input:   sequence (array of integers) to be sorted
    //          statistics describing the process of sorting the whole sequence
    // output:  sequence sorted in ascending order
    //          updated statistics
    public static void mergeSort(int[] sequence, Statistics s) {
        int length = sequence.length;
        if(length > 1) {
            int split = (int)Math.floor((double)length/2.0);
            s.memoryOps += length;
            int[] left = copyToNew(sequence, 0, split);
            int[] right = copyToNew(sequence, split, length);
            mergeSort(left, s);                                 
            mergeSort(right, s);
            merge(left, right, sequence, s);
        }
    }
    
    
    // method:  merges left and right sorted sequences into one sorted sequence.
    // input:   left and right sorted sequences
    //          statistics describing the process of sorting the whole sequence
    // output:  sorted sequence
    //          updated statistics
    public static void merge(int[] left, int [] right, int[] sequence, Statistics s) {
        int i=0, j=0, k=0;
        int leftLength = left.length;
        int rightLength = right.length;
        int sequenceLength = sequence.length;
        while(i < leftLength && j < rightLength) {
            s.comparisons++;
            s.memoryOps++;
            if(left[i] <= right[j]) {
                sequence[k] = left[i];
                i++;
            }
            else {
                sequence[k] = right[j];
                j++;
            }
            k++;
        }
        s.memoryOps += sequenceLength - k;
        if(i == leftLength)
            copyFromIdx(right, sequence, j, k);
        else
            copyFromIdx(left, sequence, i, k);
    }
    
    // method:  copies a subsequence of the source sequence into a new sequence
    // input:   source sequence (array of integers)
    //          left and right indices defining a subsequence of the source sequence
    // output:  the copy of subsequence of source sequence
    public static int[] copyToNew(int[] source, int left, int right) {
        int length = right-left;
        int[] destination = new int[length];
        for(int i = left, j=0; i<right; i++, j++)
            destination[j] = source[i];
        return destination;
    }

    // method:  copies elements from the source sequence into the destination sequence starting
    //          from indices specified for both of the sequences.
    // input:   source and destination sequences (arrays of integers)
    //          source index to start copy from
    //          destination index to start copy into
    // output:  destination sequence with some elements copied from the source sequence
    public static void copyFromIdx(int[] source, int[] destination, int srcIdx, int dstIdx) {
        int srcLen = source.length;
        int dstLen = destination.length;
        for(; srcIdx < srcLen && dstIdx < dstLen; srcIdx++, dstIdx++)
            destination[dstIdx] = source[srcIdx];
    }
    
    ////////////////////////////////// Generel Purpose Methods ///////////////////////////////////
    
    // method:  swaps two values in the sequence
    // input:   sequence (array of integers) in which the elements have to be swapped
    //          two indices fo elements which values must be swapped
    // output:  sequence with swapped elements
    public static void swap(int[] sequence, int i, int j) {
        int buffer = sequence[i];
        sequence[i] = sequence[j];
        sequence[j] = buffer;
    }
}