#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <unistd.h>

#include "utils.hpp"
#include "Trie.hpp"

using namespace std;
namespace fs = std::filesystem;

float *vb = nullptr, *vq = nullptr;
vector<vector<int>> lb, lq;
std::pair<uint32_t, float>* gt = nullptr;

int main(int argc, char* argv[]) {
    int num_base = 1000000;
    int nq = 10000;
    int dim = 128;
    int m = 32;
    int ef = 100;
    int label_num = 52;
    double alpha = 1;
    double build_ratio = 0.8;
    std::string mode = "search";
    std::string dataset = "sift";

    std::string base_path;
    std::string query_path;
    std::string base_label_path;
    std::string query_label_path;
    std::string gt_path;
    std::string query_type;
    std::string result_file;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--num_base" && i + 1 < argc) num_base = atoi(argv[++i]);
        else if (arg == "--nq" && i + 1 < argc) nq = atoi(argv[++i]);
        else if (arg == "--dim" && i + 1 < argc) dim = atoi(argv[++i]);
        else if (arg == "--m" && i + 1 < argc) m = atoi(argv[++i]);
        else if (arg == "--ef" && i + 1 < argc) ef = atoi(argv[++i]);
        else if (arg == "--label_num" && i + 1 < argc) label_num = atoi(argv[++i]);
        else if (arg == "--alpha" && i + 1 < argc) alpha = atof(argv[++i]);
        else if (arg == "--build_ratio" && i + 1 < argc) build_ratio = atof(argv[++i]);
        else if (arg == "--mode" && i + 1 < argc) mode = argv[++i];
        else if (arg == "--dataset" && i + 1 < argc) dataset = argv[++i];
        else if (arg == "--base_path" && i + 1 < argc) base_path = argv[++i];
        else if (arg == "--query_path" && i + 1 < argc) query_path = argv[++i];
        else if (arg == "--base_label_path" && i + 1 < argc) base_label_path = argv[++i];
        else if (arg == "--query_label_path" && i + 1 < argc) query_label_path = argv[++i];
        else if (arg == "--gt_path" && i + 1 < argc) gt_path = argv[++i];
        else if (arg == "--query_type" && i + 1 < argc) query_type = argv[++i];
        else if (arg == "--result_file" && i + 1 < argc) result_file = argv[++i];
    }

    if (mode != "build" && mode != "search" && mode != "insert") {
        std::cerr << "Invalid mode! Use --mode build, search, or insert" << std::endl;
        return 1;
    }

    std::string index_dir = "./index/" + dataset + "/";
    fs::create_directories(index_dir);
    std::cout << index_dir << std::endl;
    std::string index_filename = index_dir + dataset +
        "_nb" + std::to_string(num_base) +
        "_dim" + std::to_string(dim) +
        "_m" + std::to_string(m) +
        "_ef" + std::to_string(ef) +
        "_label" + std::to_string(label_num) +
        "_alpha" + std::to_string(alpha) + ".ind2";

    int build_percent = static_cast<int>(build_ratio * 100 + 0.5);
    std::string base_index_filename = index_dir + dataset +
        "_nb" + std::to_string(num_base) +
        "_dim" + std::to_string(dim) +
        "_m" + std::to_string(m) +
        "_ef" + std::to_string(ef) +
        "_label" + std::to_string(label_num) +
        "_alpha" + std::to_string(alpha) +
        "_build" + std::to_string(build_percent) + ".ind2";

    if (mode == "build") {
        if (base_path.empty() || query_path.empty() || base_label_path.empty() || query_label_path.empty()) {
            std::cerr << "Please provide --base_path --query_path --base_label_path --query_label_path" << std::endl;
            return 1;
        }
        load_data(base_path.c_str(), vb, num_base, dim);
        load_data(query_path.c_str(), vq, nq, dim);
        load_label(base_label_path.c_str(), lb, num_base);
        load_label(query_label_path.c_str(), lq, nq);

        int *keylist = new int[num_base];
        for (int i = 0; i < num_base; i++) keylist[i] = i;

        TrieIndex trieIndex(dim, num_base, num_base, vb, keylist, lb, m, ef, label_num, alpha);
        trieIndex.saveIndex(index_filename);

        std::cout << "Index built and saved to: " << index_filename << std::endl;
        delete[] keylist;
        delete[] vb;
        delete[] vq;
        return 0;
    }

    auto run_search = [&](TrieIndex* trieIndex,
                          const std::string& result_file_path,
                          const std::string& gt_path_value,
                          bool include_insert_metrics,
                          double insert_time,
                          double insert_qps,
                          int build_percent_value,
                          size_t insert_count) -> int {
        std::ofstream ofs(result_file_path, std::ios::app);
        if (!ofs.is_open()) {
            std::cerr << "Failed to open result file: " << result_file_path << std::endl;
            return 1;
        }

        vector<int> efs = {10,15,20,40,60,80,100,200,400,800,1600,3200,6400,12800,25600,51200};

        vector<int> visit(num_base, -1);
        for(int efv : efs){
            double recall = 0, comp = 0, time = 0, time2 = 0;
            double avg_size = 0;
            double avg_group_size = 0;
            double avg_part_size = 0;
            double avg_part_count = 0;
            double avg_split_size = 0;
            for(int i = 0; i < num_base; i++) visit[i] = -1;
            for(int i = 0; i < nq; i++){
                vector<TrieNode*> nodeset;
                auto start = std::chrono::high_resolution_clock::now();
                if(query_type == "containment"){
                    trieIndex->superset_search(lq[i],nodeset);
                    avg_group_size += nodeset.size();
                    double size_sum = 0;
                    for(auto node: nodeset) size_sum += node->size;
                    avg_size += size_sum/num_base;

                    sort(nodeset.begin(),nodeset.end(),
                        [](const TrieNode* a, const TrieNode* b) { return a->size > b->size; });
                    int part_count = 0;
                    double part_size = 0;
                    vector<TrieNode*> nodeset1;
                    for(auto node: nodeset) {
                        part_size += node->size;
                        part_count++;
                        nodeset1.push_back(node);
                        if(part_size/size_sum >= 0.8) break;
                    }


                    trieIndex->tag ++;
                    part_count =  min(part_count,min(efv,(int)nodeset1.size()));
                    nodeset1.resize(part_count);
                    auto result = trieIndex->searchKnn(vq+i*dim, nodeset1, 10, efv, lq[i],query_type);

                    vector<TrieNode*> nodeset2;
                    double re_size = size_sum-part_size;
                    for(int i = nodeset1.size(); i < nodeset.size(); i++){
                        if(nodeset[i]->visit!=trieIndex->tag)
                            nodeset2.push_back(nodeset[i]);
                        if(nodeset2.size()>=efv) break;
                    }

                    avg_part_size += part_size/size_sum;
                    avg_part_count += part_count;
                    if(part_count<nodeset.size())avg_split_size += nodeset[part_count]->size/size_sum;

                    auto result2 = trieIndex->searchKnn(vq+i*dim, nodeset2, 10, efv, lq[i],query_type);

                    auto end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = end - start;
                    time += elapsed.count();


                    comp += trieIndex->comp;
                    time2 += trieIndex->time;
                    int x = 0;

                    while(!result.empty()){
                        int id = result.top().second;
                        visit[id]=i;
                        result.pop();
                        if(++x==10) break;
                    }
                    while(!result2.empty()){
                        int id = result2.top().second;
                        visit[id]=i;
                        result2.pop();
                        if(++x==10) break;
                    }
                }
                else if (query_type == "overlap")
                {
                    trieIndex->overlap_search(lq[i],nodeset);
                    avg_group_size += nodeset.size();
                    double size_sum = 0;
                    for(auto node: nodeset) size_sum += node->size;
                    avg_size += size_sum/num_base;

                    sort(nodeset.begin(),nodeset.end(),
                        [](const TrieNode* a, const TrieNode* b) { return a->size > b->size; });
                    int part_count = 0;
                    double part_size = 0;

                    vector<TrieNode*> nodeset1;
                    for(auto node: nodeset) {
                        part_size += node->size;
                        part_count++;
                        nodeset1.push_back(node);
                        if(part_size/size_sum >= 0.8) break;
                    }


                    trieIndex->tag ++;
                    part_count =  min(part_count,min(efv,(int)nodeset1.size()));
                    nodeset1.resize(part_count);
                    auto result = trieIndex->searchKnn(vq+i*dim, nodeset1, 10, efv, lq[i], query_type);
                    vector<TrieNode*> nodeset2;
                    double re_size = size_sum-part_size;
                    for(int i = nodeset1.size(); i < nodeset.size(); i++){
                        if(nodeset[i]->visit!=trieIndex->tag)
                            nodeset2.push_back(nodeset[i]);
                        if(nodeset2.size()>=efv) break;
                    }

                    avg_part_size += part_size/size_sum;
                    avg_part_count += part_count;
                    if(part_count<nodeset.size())avg_split_size += nodeset[part_count]->size/size_sum;

                    auto result2 = trieIndex->searchKnn(vq+i*dim, nodeset2, 10, efv, lq[i], query_type);

                    auto end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = end - start;
                    time += elapsed.count();


                    comp += trieIndex->comp;
                    time2 += trieIndex->time;
                    int x = 0;

                    while(!result.empty()){
                        int id = result.top().second;
                        visit[id]=i;
                        result.pop();
                        if(++x==10) break;
                    }
                    x = 0;
                    while(!result2.empty()){
                        int id = result2.top().second;
                        visit[id]=i;
                        result2.pop();
                        if(++x==10) break;
                    }
                }
                else if(query_type =="equality"){
                    trieIndex->equal_search(lq[i],nodeset);
                    trieIndex->tag ++;

                    auto result = trieIndex->searchKnn(vq+i*dim, nodeset, 10, efv, lq[i], query_type);
                    auto end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = end - start;
                    time += elapsed.count();


                    comp += trieIndex->comp;
                    time2 += trieIndex->time;
                    int x = 0;

                    while(!result.empty()){
                        int id = result.top().second;
                        visit[id]=i;
                        result.pop();
                        if(++x==10) break;
                    }

                }
                double recall_tmp = 0;
                for(int j = 0; j < 10; j++)
                    if(visit[gt[i*10+j].first] == i)
                        recall += 0.1/nq, recall_tmp += 0.1;
            }
            if (include_insert_metrics) {
                ofs << gt_path_value << "," << query_type << "," << build_percent_value
                    << "," << insert_count << "," << insert_time << "," << insert_qps
                    << "," << efv << "," << nq/time << "," << recall << std::endl;
                std::cout << gt_path_value << "," << query_type << "," << build_percent_value
                    << "," << insert_count << "," << insert_time << "," << insert_qps
                    << "," << efv << "," << nq/time << "," << recall << std::endl;
            } else {
                ofs << gt_path_value << "," << query_type << "," << efv << "," << nq/time << "," << recall << std::endl;
                std::cout << gt_path_value << "," << query_type << "," << efv << "," << nq/time << "," << recall << std::endl;
            }
            if(recall > 0.995 || nq/time < 1) break;
        }
        ofs.close();
        std::cout << "Results saved to: " << result_file_path << std::endl;
        return 0;
    };

    if (mode == "insert") {
        if (base_path.empty() || query_path.empty() || base_label_path.empty() || query_label_path.empty() || gt_path.empty()) {
            std::cerr << "Please provide --base_path --query_path --base_label_path --query_label_path --gt_path for insert mode" << std::endl;
            return 1;
        }
        if (build_ratio <= 0.0 || build_ratio >= 1.0) {
            std::cerr << "Invalid build ratio! Use --build_ratio between 0 and 1 (e.g., 0.8 or 0.9)" << std::endl;
            return 1;
        }

        load_data(base_path.c_str(), vb, num_base, dim);
        load_data(query_path.c_str(), vq, nq, dim);
        load_label(base_label_path.c_str(), lb, num_base);
        load_label(query_label_path.c_str(), lq, nq);

        gt = new std::pair<uint32_t, float>[nq * 10];
        load_gt_file(gt_path.c_str(), gt, nq, 10);

        size_t build_count = static_cast<size_t>(num_base * build_ratio);
        if (build_count == 0 || build_count >= static_cast<size_t>(num_base)) {
            std::cerr << "Build ratio yields invalid count. num_base=" << num_base
                      << " build_count=" << build_count << std::endl;
            delete[] vb;
            delete[] vq;
            delete[] gt;
            return 1;
        }

        TrieIndex* trieIndex = nullptr;
        bool loaded_from_disk = false;
        std::cout << "Checking for existing base index at: " << base_index_filename << std::endl;
        if (fs::exists(base_index_filename)) {
            trieIndex = TrieIndex::loadIndex(base_index_filename);
            trieIndex->alpha = alpha;
            loaded_from_disk = true;
            std::cout << "Base index loaded from disk: " << base_index_filename << std::endl;
        } else {
            int *keylist = new int[num_base];
            for (int i = 0; i < num_base; i++) keylist[i] = i;
            auto lb_s = lb;
            lb_s.resize(build_count);
            trieIndex = new TrieIndex(dim, build_count, num_base, vb, keylist, lb_s, m, ef, label_num, alpha);
            trieIndex->saveIndex(base_index_filename);
            std::cout << "Base index built and saved to: " << base_index_filename << std::endl;
            delete[] keylist;
            std::cout << "Restarting to perform insert with the saved base index." << std::endl;
            std::vector<char*> exec_args;
            exec_args.reserve(static_cast<size_t>(argc) + 1);
            for (int i = 0; i < argc; ++i) {
                exec_args.push_back(argv[i]);
            }
            exec_args.push_back(nullptr);
            execv(exec_args[0], exec_args.data());
            std::perror("execv");
            return 1;
        }

        trieIndex->alpha = alpha;
        trieIndex->max_label = trieIndex->label_num;

        if (loaded_from_disk) {
            if (trieIndex->maxNum != static_cast<size_t>(num_base)) {
                std::cerr << "Base index maxNum mismatch. expected=" << num_base
                          << " got=" << trieIndex->maxNum << std::endl;
                delete[] vb;
                delete[] vq;
                delete[] gt;
                delete trieIndex;
                return 1;
            }
            size_t old_count = trieIndex->eleCount;
            size_t vec_bytes = trieIndex->maxNum * trieIndex->data_size_;
            char* new_vec = new char[vec_bytes];
            memcpy(new_vec, trieIndex->vecData_, old_count * trieIndex->data_size_);
            delete[] trieIndex->vecData_;
            trieIndex->vecData_ = new_vec;

            int* new_key = new int[trieIndex->maxNum];
            memcpy(new_key, trieIndex->keyList_, old_count * sizeof(int));
            delete[] trieIndex->keyList_;
            trieIndex->keyList_ = new_key;
        }

        size_t insert_count = static_cast<size_t>(num_base) - trieIndex->eleCount;
        auto insert_start = std::chrono::high_resolution_clock::now();
        for (size_t i = trieIndex->eleCount; i < static_cast<size_t>(num_base); i++) {
            trieIndex->keyList_[i] = static_cast<int>(i);
            trieIndex->key2Id[static_cast<int>(i)] = static_cast<int>(i);
            char* data_ptr = reinterpret_cast<char*>(vb + i * dim);
            trieIndex->insertElement(static_cast<int>(i), lb[i], data_ptr);
        }
        auto insert_end = std::chrono::high_resolution_clock::now();
        trieIndex->trie_label_load(trieIndex->root);
        std::chrono::duration<double> insert_elapsed = insert_end - insert_start;
        double insert_time = insert_elapsed.count();
        double insert_qps = insert_time > 0 ? insert_count / insert_time : 0.0;
        std::cout << "Insert done. build_ratio=" << build_ratio
                  << " insert_count=" << insert_count
                  << " insert_time=" << insert_time
                  << " insert_qps=" << insert_qps << std::endl;

        int search_status = run_search(trieIndex, result_file, gt_path, true,
                                       insert_time, insert_qps, build_percent, insert_count);

        delete[] vb;
        delete[] vq;
        delete[] gt;
        delete trieIndex;
        return search_status;
    }

    if (query_path.empty() || query_label_path.empty() || gt_path.empty()) {
        std::cerr << "Please provide --query_path --query_label_path --gt_path for search mode" << std::endl;
        return 1;
    }
    load_data(query_path.c_str(), vq, nq, dim);
    load_label(query_label_path.c_str(), lq, nq);

    gt = new std::pair<uint32_t, float>[nq * 10];
    load_gt_file(gt_path.c_str(), gt, nq, 10);
    std::cout << index_filename << std::endl;
    auto preMem = get_memory_usage();
    TrieIndex* trieIndex = TrieIndex::loadIndex(index_filename);
    auto postMem = get_memory_usage();
    print_memory_diff("Index loading", preMem, postMem);

    int search_status = run_search(trieIndex, result_file, gt_path, false, 0.0, 0.0, 0, 0);
    delete[] vq;
    delete[] gt;
    delete trieIndex;
    return search_status;
}
