import numba
import numpy as np

# nopython=True means an error will be raised
# if fast compilation is not possible.
@numba.jit(nopython=True)
def fill_counts(row_centers,col_centers,row_indices,col_indices):
    """
          given bincenters for each row and column, and row_indices and
          col_indices from searchsorted that give the row,col index of
          each binned data point, return the counts in each 2d bin_col

          input:  row_centers, col_centers   -- vectors giving the location
                  of the center of each bin_col

                  row_indices, col_indices  -- vectors giving the bin that
                  every data point belongs in

                  hist_array  -- the output array holding the counts -- this
                   needs to be declared outside

          output:  hist_array  -- 2d array of shape [len(row_centers),len(col_centers)]
                  containing the number of datapoints that are in each row,column bin_column
                  If there are no datapoints in a bin then the bin contains np.nan
    """
    num_colbins=int(col_centers.shape[0])
    num_rowbins=int(row_centers.shape[0])
    hist_array=np.zeros((num_rowbins,num_colbins))
    num_y=row_indices.shape[0]
    for n in range(num_y): #row_indices and col_indices both size of raw data
        if col_indices[n] > 0 and row_indices[n] > 0 and \
            col_indices[n] <= num_colbins and row_indices[n] <= num_rowbins:
            bin_row=row_indices[n]-1 # '-1' to get the index of the bin center
            bin_col=col_indices[n]-1
            hist_array[bin_row, bin_col] += 1
    rows,cols=hist_array.shape
    for row in range(rows):
        for col in range(cols):
            if hist_array[row,col] < 1.:
                hist_array[row,col]=np.nan
    return hist_array

            
def hist2d(col_raw,row_raw,col_edges,row_edges):
    """
      Produce a 2-d histogram (for example, of temperature (y) vs.
        vertical velocity (x) data)  binned into temperature,wvel
        bins
      input: row_raw,col_raw: data vectors of the row variable (temperature)
             and the column variable (wvel)
             col_edges, row_edges:  coordinates of the bin edges for each variables
      returns:  counts,col_centers,row_centers

      Example, given 10,000 temperature measurements to be binned into 20 bins, and
               20,000 wvel measurements to be binned into 10 bins, return
               counts as a [20,10]  array with the number of measurements that fall
               into each bin
    """
    col_centers=(col_edges[:-1] + col_edges[1:])/2.
    row_centers=(row_edges[:-1] + row_edges[1:])/2.
    col_indices=np.asarray(np.searchsorted(col_edges, col_raw.flat, 'right'),dtype=np.int64)
    row_indices=np.asarray(np.searchsorted(row_edges, row_raw.flat, 'right'),dtype=np.int64)
    hist_array=fill_counts(row_centers,col_centers,row_indices,col_indices)
    return hist_array,col_centers,row_centers

