/**
 * @author      Bogdan Kovch
 * @class       CS 350
 * @assignment  final project
 * @file        Statistics.java
 * @date        May 29, 2013 - June 6, 2013
 */

/**
 * This file defines the class used to populate instances encapsulating statistics data about
 * sorting algorithm performance.
 */
    
package sortalgs;

public class Statistics {
    public int length = 0;          // number of elements in the sequence
    public long comparisons = 0;    // number of element comparisons an algorithm performed
    public long memoryOps = 0;      // number of memory operations; that is the number of
                                    // writes into sequence elements 
    public long time = 0;           // total time taken to sort the sequence in nanoseconds
    
    Statistics(int length){         // constructor defines the length of the sequence
        this.length = length;
    }
}
