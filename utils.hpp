#ifndef LABELFILTEREDANNS_UTILS_HPP
#define LABELFILTEREDANNS_UTILS_HPP
#include <omp.h>
#include <fstream>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <chrono>
#include <algorithm>

using namespace std;

void load_data(const char* filename, float*& data, int num, int dim) {
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cout << "open file error" << std::endl;
        exit(-1);
    }
    in.read((char*)&dim, 4);
    in.seekg(0, std::ios::end);
    std::ios::pos_type ss = in.tellg();
    size_t fsize = (size_t)ss;
    num = (unsigned)(fsize / (dim + 1) / 4);
    data = new float[(size_t)num * (size_t)dim];

    in.seekg(0, std::ios::beg);
    for (size_t i = 0; i < num; i++) {
        in.seekg(4, std::ios::cur);
        in.read((char*)(data + i * dim), dim * 4);
    }
    in.close();
}

void load_label(const char* filename, vector<vector<int> >& data, int num) {
    std::string line;
    size_t lines_read = 0;

    std::ifstream file(filename);

    data.reserve(num > 0 ? num : 100);
    while ((num == 0 || lines_read < num) &&
           std::getline(file, line))
    {
        if (line.empty()) continue;

        std::vector<int> row;
        std::stringstream ss(line);
        std::string cell;

        size_t comma_count = std::count(line.begin(), line.end(), ',');
        row.reserve(comma_count + 1);

        while (std::getline(ss, cell, ',')) {
            try {
                row.push_back(std::stoi(cell));
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "Error: Invalid number format at line " +
                    std::to_string(lines_read + 1) +
                    ", value: '" + cell + "'"
                );
            }
        }

        data.push_back(std::move(row));
        lines_read++;
    }
}

float getDistance(float* a, float *b, int dim){
    float ans = 0;
    for(int i = 0; i < dim; i++)
        ans+=(a[i] - b[i]) * (a[i] - b[i]);
    return ans;
}

bool check(vector<int> a, vector<int>b){
    for(auto x:a){
        bool flag = false;
        for(auto y:b){
            if(x==y) flag = true;
        }
        if(!flag)return false;
    }
    return true;
}

void load_gt_file(const std::string& filename, std::pair<uint32_t, float>* gt, uint32_t num_queries, uint32_t K) {
        std::ifstream fin(filename, std::ios::binary);
        fin.read(reinterpret_cast<char*>(gt), num_queries * K * sizeof(std::pair<uint32_t, float>));
        std::cout << "Ground truth loaded from " << filename << std::endl;
    }

void get_gt(float* vb, float* vq, vector<vector<int> > &lb, vector<vector<int> >& lq, int numb, int numq, int k,int dim, ofstream &outfile){
    for(int id = 0; id < numq; id++){
        vector<pair<float,int> > result;
        for(int i = 0; i < numb; i++){
            if(check(lq[id],lb[i]))
                result.push_back({getDistance(vb+i*dim,vq+id*dim, dim),i});
        }
        sort(result.begin(), result.end());
        result.resize(k);
        for(int i = 0; i < k; i++)
            outfile<<result[i].second<<' ';
        outfile<<endl;
    }
}

size_t get_memory_usage() {
    std::ifstream status_file("/proc/self/status");
    std::string line;
    size_t memory_usage = 0;

    while (std::getline(status_file, line)) {
        if (line.find("VmRSS:") != std::string::npos) {
            size_t pos = line.find_first_of("0123456789");
            if (pos != std::string::npos) {
                memory_usage = std::stoul(line.substr(pos));
            }
            break;
        }
    }
    return memory_usage;
}

void print_memory_diff(const std::string& label,
                       size_t before,
                       size_t after) {
    std::cout << "\n----- " << label << " -----\n";
    std::cout << "Memory before: " << before << " KB\n";
    std::cout << "Memory after: " << after << " KB\n";
    std::cout << "Diff: " << (float)(static_cast<long>(after) - static_cast<long>(before))/1024
              << " MB\n";
    std::cout << "--------------------------\n";
}


#endif //LABELFILTEREDANNS_UTILS_HPP
