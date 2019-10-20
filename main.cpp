/* motif-search project for CSCI 590: introduction to Bioinfomatics
 * Connor Horne
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <filesystem>
#include <unistd.h>

namespace fs = std::filesystem;

// score a position and return the score and consensus motif
std::pair<std::string, int> score(std::vector<int> &position, int index,
                                  std::vector<std::string> &sequences,
                                  int width) {

  // Matrix, the format is as follows
  // 0 = A
  // 1 = C
  // 2 = G
  // 3 = T
  // static because it was getting allocated on each run which didn't leak but
  // was expensive so static
  static std::vector<std::vector<int>> consensus_matrix(width,
                                                        std::vector<int>(4, 0));

  // initialize it out on a new run
  for (auto &col : consensus_matrix) {
    std::fill(col.begin(), col.end(), 0);
  }

  // loop through sequences
  for (int i = 0; i < index; i++) {

    // get postion
    int location = position.at(i);

    // loop through possible motif
    for (int j = 0; j < width; j++) {
      int motif_index = j + location;
      auto amino = sequences.at(i).at(motif_index);

      switch (amino) {
      case 'a':
        consensus_matrix.at(j).at(0) = consensus_matrix.at(j).at(0) + 1;
        break;
      case 'c':
        consensus_matrix.at(j).at(1) = consensus_matrix.at(j).at(1) + 1;
        break;
      case 'g':
        consensus_matrix.at(j).at(2) = consensus_matrix.at(j).at(2) + 1;
        break;
      case 't':
        consensus_matrix.at(j).at(3) = consensus_matrix.at(j).at(3) + 1;
        break;
      }
    }
  }

  int score = 0;
  std::string motif = "";

  // get the total score & consensus motif from the consensus matrix
  for (const auto &col : consensus_matrix) {
    auto result_iter = std::max_element(col.begin(), col.end());
    auto result = std::distance(col.begin(), result_iter);

    // create motif from consensus matix
    switch (result) {
    case 0:
      motif.push_back('a');
      break;
    case 1:
      motif.push_back('c');
      break;
    case 2:
      motif.push_back('g');
      break;
    case 3:
      motif.push_back('t');
    }

    score = score + col.at(result);
  }

  return std::make_pair(motif, score);
}

// jump down level, traverse sublevel, or jump to upperlevel
int next_vertex(std::vector<int> &position, int index, int length,
                int k_lmers) {

  // jump to lower level
  if (index < length) {

    position.at(index + 1) = 1;

    return index + 1;
  } else {

    for (int j = length; j >= 0; j--) {

      if (position.at(j) < k_lmers) {
        position.at(j) = position.at(j) + 1;
        return j;
      }
    }
  }

  return -1;
}

// jump past a branch if it sucks
int bypass(std::vector<int> &position, int index, int k_lmers) {

  for (int j = index; j >= 0; j--) {

    if (position.at(j) < k_lmers) {
      position.at(j) = position.at(j) + 1;
      return j;
    }
  }

  return -1;
}

// branch and bound algorithm from intro. to bioinformatics
// Jones & Pevzner pg. 111
std::pair<std::string, int> branch_bound(std::vector<std::string> &sequences,
                                         int width) {

  std::vector<int> position(sequences.size(), 0);

  int length = sequences.at(0).size();
  int num_seqs = sequences.size() - 1;
  int index = 0;
  int k_lmers = length - width;

  int best_score = 0;
  std::string best_motif = "";

  while (index > -1) {
    if (index < num_seqs) {

      auto score_result = score(position, index, sequences, width);

      int optimistic_score = score_result.second + ((num_seqs - index) * width);

      if (optimistic_score <= best_score) {
        index = bypass(position, index, k_lmers);
      } else {
        index = next_vertex(position, index, num_seqs, k_lmers);
      }
    } else {
      auto score_result = score(position, index, sequences, width);

      if (score_result.second > best_score) {
        best_score = score_result.second;
        best_motif = score_result.first;
      }
      index = next_vertex(position, index, num_seqs, k_lmers);
    }
  }

  return std::make_pair(best_motif, best_score);
}

// modified greedy approach to motif search from intro. to bioinfomatics
// Jones & Pevzner pg. 136
std::pair<std::string, double> greedy(std::vector<std::string> &sequences,
                                      int width) {

  int length = sequences.at(0).size();
  int num_seqs = sequences.size();
  int k_lmers = length - width;

  // get seed sequences
  int seed_1 = 0;
  int seed_2 = 1;

  std::vector<std::pair<int, int>> rand_seq_pairs;

  // get random setup
  std::random_device rd;

  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dist(0, num_seqs - 1);

  int swap_1 = -1;
  int swap_2 = -1;

  // create 20 pairs of random sequences
  for (int i = 0; i < 20; i++) {

    swap_1 = -1;
    swap_2 = -1;

    // make sure not to get a same swap
    while ((swap_1 == -1) || (swap_2 == -1)) {
      swap_1 = dist(gen);
      swap_2 = dist(gen);

      if (swap_1 == swap_2) {
        swap_1 = -1;
        swap_2 = -1;
      }
    }
    rand_seq_pairs.push_back(std::make_pair(swap_1, swap_2));
  }

  std::vector<int> position(sequences.size(), 0);
  std::vector<int> best_motif(sequences.size(), 0);

  std::vector<std::string> rand_sequences(2);
  std::pair<int, int> best_pair(0, 0);

  int best_score = 0;

  // loop through the random generated sequences
  // and check for best pair
  for (const auto &pair : rand_seq_pairs) {

    // gah copying too much
    rand_sequences.at(0) = sequences.at(pair.first);
    rand_sequences.at(1) = sequences.at(pair.second);

    std::vector<int> rand_position(2, 0);
    std::vector<int> rand_best_motif(2, 0);

    // create seed matrix
    for (int i = 0; i < k_lmers; i++) {
      rand_position.at(seed_1) = i;
      for (int j = 0; j < k_lmers; j++) {
        rand_position.at(seed_2) = j;

        auto position_result = score(rand_position, 2, rand_sequences, width);
        auto best_motif_result =
            score(rand_best_motif, 2, rand_sequences, width);

        // update if this is the best lmer this pair of sequences
        if (position_result.second > best_motif_result.second) {
          rand_best_motif.at(seed_1) = rand_position.at(seed_1);
          rand_best_motif.at(seed_2) = rand_position.at(seed_2);

          // update if this is the best l-mer we've seen
          if (position_result.second > best_score) {
            best_pair = pair;
            best_score = position_result.second;

            best_motif.at(seed_1) = rand_best_motif.at(seed_1);
            best_motif.at(seed_2) = rand_best_motif.at(seed_2);
          }
        }
      }
    }
  }

  // swap for running the algorithm correctly as I am too lazy to write multiple
  // scoring functions
  sequences.at(0).swap(sequences.at(best_pair.first));
  sequences.at(1).swap(sequences.at(best_pair.second));

  position.at(seed_1) = best_motif.at(seed_1);
  position.at(seed_2) = best_motif.at(seed_2);

  // extend seed matrix to all sequences
  for (int i = 2; i < num_seqs; i++) {
    for (int j = 0; j < k_lmers + 1; j++) {

      position.at(i) = j;

      auto position_result = score(position, i + 1, sequences, width);
      auto best_motif_result = score(best_motif, i + 1, sequences, width);

      if (position_result.second > best_motif_result.second) {
        best_motif.at(i) = position.at(i);
      }
    }
    position.at(i) = best_motif.at(i);
  }

  auto final_motif = score(position, num_seqs, sequences, width);

  return std::make_pair(final_motif.first, final_motif.second);
}

void help() {
  std::cout << "MotifSearch sequence_filename.txt w op \n"
            << "Description: Program for finding motifs in sequences \n"
            << "using branch and bound method and greedy search as \n"
            << "specified by project details. \n"
            << "sequence_file: file containing sequences \n"
            << "w: length of desired motif\n"
            << "op: selected method, bb for branch and bound, greedy\n"
            << "for greedy search" << std::endl;
}

// returns the sequences from a file
void parse_file(fs::path file_path, std::vector<std::string> &sequences) {

  std::ifstream file(file_path.string());

  if (file.is_open()) {
    std::string line;

    while (std::getline(file, line)) {
      // validate line?
      if (line.at(0) != '>') {
        sequences.push_back(line);
      }
    }
  }
}

int main(int argc, char **argv) {

  int c;

  // optional arguments
  while ((c = getopt(argc, argv, "w:")) != -1) {

    // handle optional arguments
    switch (c) {
    case 'h':
    case '?':
    default:
      help();
      break;
    }
  }

  // check if their are enough args
  if ((argc - optind) != 3) {
    std::cerr << "not correct number of arguments" << std::endl;
    return -1;
  }

  std::string file_name = argv[optind];
  std::string str_width = argv[optind + 1];
  std::string operation = argv[optind + 2];

  // test path
  fs::path file_path(file_name);
  if (!fs::exists(file_path)) {
    std::cerr << "file " << file_name << " does not exists" << std::endl;
  }

  // test width and get int val
  int width;
  try {
    width = std::stoi(str_width);
  } catch (const std::exception &e) {
    std::cout << "width cannot be converted to int: " << e.what() << std::endl;
    return -1;
  }

  if (operation != "bb" && operation != "greedy") {
    std::cerr << "operation is not vaild" << std::endl;
    return -1;
  }

  // readin sequences
  std::vector<std::string> sequences;
  parse_file(file_name, sequences);

  // run the algorithm and get the sequence and score
  std::pair<std::string, double> motif_result;
  if (operation == "bb") {
    motif_result = branch_bound(sequences, width);
  } else {
    motif_result = greedy(sequences, width);
  }

  std::cout << motif_result.first << "\t" << motif_result.second << std::endl;
}
