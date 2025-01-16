#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List sumMultipleRleListsOptimized2(List rleValuesList, List rleLengthsList) {
  int num_chromosomes = as<List>(rleValuesList[0]).size();  // Number of chromosomes
  int num_lists = rleValuesList.size();  // Number of RleLists

  List resultList(num_chromosomes);  // Resulting list of summed Rle values for each chromosome

  // Loop over chromosomes
  for (int chrom = 0; chrom < num_chromosomes; chrom++) {
    // Precompute the total length to reserve memory
    int estimated_size = 0;
    std::vector<NumericVector> values(num_lists);
    std::vector<IntegerVector> lengths(num_lists);

    for (int list_idx = 0; list_idx < num_lists; list_idx++) {
      List currentValuesList = rleValuesList[list_idx];
      List currentLengthsList = rleLengthsList[list_idx];

      values[list_idx] = currentValuesList[chrom];
      lengths[list_idx] = currentLengthsList[chrom];

      // Add up lengths of runs in each list to get the total number of runs
      estimated_size += values[list_idx].size();
    }

    // Allocate raw arrays to avoid push_back overhead
    double* result_values = new double[estimated_size];
    int* result_lengths = new int[estimated_size];
    int result_idx = 0;

    // Track positions in each list
    std::vector<int> positions(num_lists, 0);
    std::vector<int> current_pos(num_lists, 0);

    // Process each run
    while (true) {
      bool all_done = true;
      double sum_value = 0;
      int min_length = INT_MAX;

      // Calculate sum and minimum length
      for (int list_idx = 0; list_idx < num_lists; list_idx++) {
        if (positions[list_idx] < values[list_idx].size()) {
          all_done = false;
          sum_value += values[list_idx][positions[list_idx]];
          int available_length = lengths[list_idx][positions[list_idx]] - current_pos[list_idx];
          if (available_length < min_length) {
            min_length = available_length;
          }
        }
      }

      if (all_done) break;

      // Store the result
      if (result_idx > 0 && result_values[result_idx - 1] == sum_value) {
        result_lengths[result_idx - 1] += min_length;
      } else {
        result_values[result_idx] = sum_value;
        result_lengths[result_idx] = min_length;
        result_idx++;
      }

      // Advance positions
      for (int list_idx = 0; list_idx < num_lists; list_idx++) {
        if (positions[list_idx] < values[list_idx].size()) {
          current_pos[list_idx] += min_length;
          if (current_pos[list_idx] == lengths[list_idx][positions[list_idx]]) {
            positions[list_idx]++;
            current_pos[list_idx] = 0;
          }
        }
      }
    }

    // Convert raw arrays back to vectors of appropriate length
    NumericVector final_values(result_values, result_values + result_idx);
    IntegerVector final_lengths(result_lengths, result_lengths + result_idx);

    // Clean up memory
    delete[] result_values;
    delete[] result_lengths;

    // Store the result for the chromosome
    resultList[chrom] = List::create(
      Named("values") = final_values,
      Named("lengths") = final_lengths
    );
  }

  return resultList;
}
