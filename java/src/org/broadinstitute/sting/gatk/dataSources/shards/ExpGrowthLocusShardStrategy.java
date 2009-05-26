package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.List;

/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 8:23:19 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 6, 2009
 * <p/>
 * Class LinearShard
 * <p/>
 * A exponential strategy
 */
public class ExpGrowthLocusShardStrategy extends LocusShardStrategy {

    // fixed size
    private long baseSize = 100000;
    private long currentExp = 0;

    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param dic the seq dictionary
     */
    ExpGrowthLocusShardStrategy(SAMSequenceDictionary dic, long startSize) {
        super(dic);
        this.baseSize = startSize;
        currentExp = 0;
    }

    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param strat the shatter to convert from
     */
    ExpGrowthLocusShardStrategy(LocusShardStrategy strat) {
        super(strat);
        this.baseSize = strat.nextShardSize();
        currentExp = 0;
    }

    /**
     * The constructor, for a genomic list, start size, and a reference dictionary
     *
     * @param dic       the reference dictionary
     * @param startSize the starting size of the shard
     * @param lst       locations to iterate from
     */
    ExpGrowthLocusShardStrategy(SAMSequenceDictionary dic, long startSize, GenomeLocSortedSet lst) {
        super(dic, lst);
        this.baseSize = startSize;
        this.currentExp = 0;
    }


    /**
     * set the next shards size
     *
     * @param size adjust the next size to this
     */
    public void adjustNextShardSize(long size) {
        baseSize = size;
        currentExp = 0;
    }

    /**
     * This is how the various shards strategies implements their approach
     *
     * @return the next shard size
     */
    protected long nextShardSize() {
        // we grow the exponentially, we just have to make sure we start at zero
        ++currentExp;
        return (long) Math.floor(Math.pow((double) baseSize, (double) currentExp));
    }

}
