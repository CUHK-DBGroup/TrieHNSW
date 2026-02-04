#ifndef LABELFILTEREDANNS_TRIE_HPP
#define LABELFILTEREDANNS_TRIE_HPP

#include <unordered_map>
#include <map>
#include <vector>
#include <functional>
#include "hnswlib/hnswlib.h"

using namespace std;
using namespace hnswlib;

typedef int LabelType;

struct CompareByFirst {
    constexpr bool operator()(std::pair<float, tableint> const& a,
                              std::pair<float, tableint> const& b) const noexcept {
        return a.first < b.first;
    }
};

typedef std::priority_queue<std::pair<float, tableint>, std::vector<std::pair<float , tableint>>, CompareByFirst> ResultHeap;

using LabelMatchFunc = std::function<bool(unsigned long long, unsigned long long)>;

struct TrieNode {
    map<int, TrieNode *> child;
    TrieNode *parent, *son;
    int size, other_size;
    int actual_size, graph_size;
    int level;
    int flag;
    int visit = 0;
    int layer;
    int entry;
    vector<int> ids;
    TrieNode* pre_index;

    TrieNode(TrieNode *p, int i) {
        parent = p;
        flag = false;
        layer = i;
        size = actual_size = graph_size = 0;
        other_size = 0;
        pre_index = this;
        son = nullptr;
        entry = -1;
    }
};



class TrieIndex {

public:

    TrieIndex(): space(128) {;
        };
    TrieIndex(
            int d,
            size_t eleNum,
            size_t maxEleNum,
            float* vecData,
            int* keyList,
            vector< vector<int > > &labels,
            int m,
            int ef_con,
            int label_num_,
            double alpha_
    ):
            M(m),ef_construction(ef_con), space(d), dim(d), linklist(maxEleNum), eleCount(eleNum), maxNum(maxEleNum), max_layer(maxEleNum),
             labelhash(maxEleNum), belong(maxEleNum), max_layer_i(maxEleNum), alpha(alpha_), label_num(label_num_), hashedlabel(label_num_ + 1)
    {

        item_labels = labels;
        for(int i = 0; i < item_labels.size(); i++){
            sort(item_labels[i].begin(), item_labels[i].end());
            item_labels[i].push_back(label_num + 1);
        }

        visited_array = new unsigned int[maxEleNum];
        memset(visited_array, 0, sizeof(unsigned int) * maxEleNum);
        std::random_device rd;
        eng = std::mt19937 (rd());

        max_label = label_num;
        eng_label = std::mt19937(37);

        data_size_ = space.get_data_size();
        fstdistfunc_ = space.get_dist_func();
        dist_func_param_ = space.get_dist_func_param();

        keyList_ = new int[maxEleNum];
        valueList_ = new int[maxEleNum];
        vecData_ = (char*)(new float[maxEleNum * dim]);
        isDeleted = new bool[maxEleNum];
        memset(isDeleted,0,maxEleNum);
        memcpy(keyList_,keyList, eleNum * sizeof(int));
        memcpy(vecData_,vecData, eleNum * dim * sizeof(float));

        int random_seed = 100;
        level_generator_.seed(random_seed);
        labelhash_generator_ = std::uniform_int_distribution<int>(1, label_num);
        for (int i = 1; i <= label_num; i++) {
            if (i <= label_the) hashedlabel[i] = i-1;
            else hashedlabel[i] = labelhash_generator_(eng_label) % (64-label_the) + label_the;
        }

        mult_ = 1 / log(1.0 * M);
        revSize_ = mult_;

        sizeLinkList = (M * sizeof(tableint) + sizeof(linklistsizeint));

        space = hnswlib::L2Space(dim);

        sortedArray.reserve(eleNum);

        for(int i = 0; i < eleNum; i++) {
            insert(i, item_labels[i]);
            labelhash[i] = setBits(item_labels[i]);
        }

        for(int i = 0; i < eleNum; i++) {
            max_layer[i] = item_labels[i].size();
        }

        for(int i = 0; i < eleNum; i++){
            key2Id[keyList[i]] = i;
            max_layer_i[i] = getRandomLevel(revSize_);
            linklist[i] = (char *) malloc( (max_layer[i] + 1) * sizeLinkList * (max_layer_i[i] + 2));
            sorted_array.push_back( {-1,i} );
        }

        build_index(root,0, eleNum - 1, 0, 0);

    }


    void insert(int id, const vector<int> &label_set) {
        TrieNode *node = root;
        if(root == nullptr) node = root = new TrieNode(nullptr, 0);
        for (int i = 0; i < label_set.size(); i++) {
            if (node->child.count(label_set[i]) == 0) {
                node->child[label_set[i]] = new TrieNode(node, i + 1);
            }
            node = node->child[label_set[i]];
        }
        node->ids.push_back(id);
    }

    void insert2(int id, const vector<int> &label_set) {
        TrieNode *node = root;
        if(root == nullptr) node = root = new TrieNode(nullptr, 0);
        for (int i = 0; i < label_set.size(); i++) {
            node->actual_size++;
            node->graph_size++;
            node->size++;
            if (node->child.count(label_set[i]) == 0) {
                node->child[label_set[i]] = new TrieNode(node, node->layer + 1);
                if(node->size == 0){
                    node->child[label_set[i]]->pre_index = node->pre_index;
                    node->son = node->child[label_set[i]];
                }
                else node->child[label_set[i]]->pre_index = node->child[label_set[i]];
            }
            if(node->child[label_set[i]]!=node->son)
                node->other_size++;
            node = node->child[label_set[i]];
        }
        node->ids.push_back(id);
        node->actual_size++;
        node->graph_size++;
        node->size++;
        max_layer[id] = node->pre_index->layer;
    }

    void insertElement(int id, const vector<int> &label_set, char* data){
        item_labels[eleCount] = label_set;
        sort(item_labels[eleCount].begin(),item_labels[eleCount].end());
        item_labels[eleCount].push_back(max_label + 1);
        labelhash[eleCount] = setBits(item_labels[eleCount]);
        memcpy(vecData_+ dim * sizeof(float) * eleCount, data, dim * sizeof(float));
        insert2(id,item_labels[eleCount]);
        belong[eleCount].resize(max_layer[eleCount] + 1);
        max_layer_i[eleCount] = getRandomLevel(revSize_);
        if(alpha==2)
            adjustTrie(id, item_labels[eleCount]);
        max_layer[eleCount] = getLayer(root, item_labels[eleCount]);
        belong[eleCount].resize( max_layer[eleCount] + 1);
        linklist[eleCount] = ((char *) malloc( (max_layer[eleCount] + 1) * sizeLinkList * (max_layer_i[eleCount] + 2)));
        insertIntoGraph(root, id, item_labels[eleCount], data,0);

        eleCount ++;
    }

    int getLayer(TrieNode *root, const vector<int> &label_set) {
        TrieNode *node = root;
        for (int i = 0; i < label_set.size(); i++) {
            node = node->child[label_set[i]];
        }
        return node->pre_index->layer;
    }

    void deleteElement(int id){
        auto label_set = item_labels[id];
        isDeleted[id] = true;
        deleteFromTrie(id, label_set);
        max_layer_i[eleCount] = getRandomLevel(revSize_);
    }

    void getbelong(int id, const vector<int> &label_set) {
        TrieNode *node = root;
        for (int i = 0; i < label_set.size(); i++) {
            if(node->pre_index == node){
                belong[id].resize(node->layer + 1);
                belong[id][node->layer] = node;
            }
            node = node->child[label_set[i]];
        }
        if(node->pre_index == node){
                belong[id].resize(node->layer + 1);
                belong[id][node->layer] = node;
            }
    }

    bool isOverlap(const vector<int> &q, const vector<int> &a) const{
        for(auto x:q){
            for(auto y:a){
                if(y > label_num) continue;
                if(x==y) return true;
            }
        }
        return false;
    }

    bool isEquality(const vector<int> &q, const vector<int> &a) const{
        int valid_size = 0;
        for(auto y : a){
            if(y <= label_num) valid_size++;
        }
        if(q.size()!=valid_size) return false;
        for(auto x:q){
            bool flag = false;
            for(auto y:a){
                if(y > label_num) continue;
                if(x==y) {
                    flag = true;
                    break;
                }
            }
            if(!flag) return false;
        }
        return true;
    }

    bool isContainment(const vector<int> &q, const vector<int> &a) const{
        for(auto x:q){
            bool flag = false;
            for(auto y:a){
                if(x==y) {
                    flag = true;
                    break;
                }
            }
            if(!flag) return false;
        }
        return true;
    }


    std::priority_queue<std::pair<float, hnswlib::labeltype>> searchKnn(float *vecData, vector<TrieNode *> &node_set, int k,int ef_s, vector<int> &x, string type){
        vector<tableint> ep_ids;
        std::priority_queue<std::pair<float, hnswlib::labeltype>> top;
        unsigned long long ql = setBits(x);
        comp = 0;
        for (auto node:node_set){
            if(node->visit == tag) continue;
            ep_ids.push_back(findEntry(vecData, node, node->entry));
        }
        ResultHeap result;
        if(type == "containment") result = searchBaseLayer0(ep_ids,vecData,ql,x,ef_s);
        else if (type == "equality") result = searchBaseLayer0forEquality(ep_ids,vecData,ql,x,ef_s);
        else if(type == "overlap") result = searchBaseLayer0ForOverlap(ep_ids,vecData,ql,x,ef_s);
        while(result.size()) {
            auto [dist, id] = result.top();
            result.pop();
            if (label_num > label_the) {
                if (type == "containment" && !isContainment(x, item_labels[id])) continue;
                else if (type == "overlap" && !isOverlap(x, item_labels[id])) continue;
                else if (type == "equality" && !isEquality(x, item_labels[id])) continue;
            }
            top.push({dist, keyList_[id]});
        }
        while(top.size() > k) top.pop();
        return top;
    }

    std::priority_queue<std::pair<float, hnswlib::labeltype>> searchKnn2(float *vecData, vector<TrieNode *> &node_set, int k,int ef_s, vector<int> &x, string type){
        vector<tableint> ep_ids;
        std::priority_queue<std::pair<float, hnswlib::labeltype>> top;
        unsigned long long ql = setBits(x);
        comp = 0;
        time = 0;
        for (auto node:node_set){
            ep_ids.clear();
            ep_ids.push_back(node->entry);

            ResultHeap result = searchBaseLayer0(ep_ids,vecData,ql,x,ef_s, node->layer);

            while(result.size() > k) result.pop();

            while(!result.empty()){
                auto r = result.top();
                result.pop();
                top.push({-r.first,keyList_[r.second]});
            }
        }
        return top;
    }


    void trie_label_load(TrieNode *node) {
        for (auto ptr: node->child) {
            if (ptr.first < 0 || ptr.first >= label_num) break;
            label_node[ptr.first].push_back(ptr.second);
            trie_label_load(ptr.second);
        }
    }

    void superset_search(const vector<int> &label_set,
                         vector<TrieNode *> &node_set) {
        tag2++;
        superset_search(root, label_set, 0, node_set);
    }

    int superset_search(const vector<int> &label_set) {
        return superset_search(root, label_set, 0);
    }

    void subset_search(const vector<int> &label_set,
                         vector<TrieNode *> &node_set) {
        subset_search(root, label_set, 0, node_set);
    }

    void equal_search(const vector<int> &label_set,
                       vector<TrieNode *> &node_set) {
        tag2++;
        equal_search(root, label_set, 0, node_set);
    }

    void overlap_search(const vector<int> &label_set,
                      vector<TrieNode *> &node_set) {
        tag2++;
        overlap_search(root, label_set, 0, node_set);
    }


    void saveIndex(const string& filename) {
        ofstream out(filename, ios::binary);
        if (!out) {
            throw runtime_error("Cannot open file for writing: " + filename);
        }

        // 1. Save basic parameters
        out.write(reinterpret_cast<const char*>(&M), sizeof(M));
        out.write(reinterpret_cast<const char*>(&ef_construction), sizeof(ef_construction));
        out.write(reinterpret_cast<const char*>(&dim), sizeof(dim));
        out.write(reinterpret_cast<const char*>(&maxNum), sizeof(maxNum));
        out.write(reinterpret_cast<const char*>(&eleCount), sizeof(eleCount));
        out.write(reinterpret_cast<const char*>(&mult_), sizeof(mult_));
        out.write(reinterpret_cast<const char*>(&revSize_), sizeof(revSize_));
        out.write(reinterpret_cast<const char*>(&label_num), sizeof(label_num));

        // 2. Save raw data arrays
        out.write(vecData_, eleCount * dim * sizeof(float));
        out.write(reinterpret_cast<const char*>(keyList_), eleCount * sizeof(int));
        out.write(reinterpret_cast<const char*>(isDeleted), maxNum * sizeof(bool));

        // 3. Save item labels
        size_t num_items = item_labels.size();
        out.write(reinterpret_cast<const char*>(&num_items), sizeof(num_items));
        for (const auto& labels : item_labels) {
            size_t size = labels.size();
            out.write(reinterpret_cast<const char*>(&size), sizeof(size));
            out.write(reinterpret_cast<const char*>(labels.data()), size * sizeof(int));
        }

        // max_layer, max_layer_i, labelhash
        out.write(reinterpret_cast<const char*>(max_layer.data()), maxNum * sizeof(int));
        out.write(reinterpret_cast<const char*>(max_layer_i.data()), maxNum * sizeof(int));
        out.write(reinterpret_cast<const char*>(labelhash.data()), maxNum * sizeof(unsigned long long));

        // 4. Save Trie tree structure
        saveTrieNode(out, root);

        // 5. Save key2Id mapping
        size_t map_size = key2Id.size();
        out.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
        for (const auto& pair : key2Id) {
            out.write(reinterpret_cast<const char*>(&pair.first), sizeof(int));
            out.write(reinterpret_cast<const char*>(&pair.second), sizeof(int));
        }

        // 6. Save link lists
        for (size_t id = 0; id < eleCount; id++) {
            int max_layer_val = max_layer[id];
            int max_layer_i_val = max_layer_i[id];
            out.write(reinterpret_cast<const char*>(&id), sizeof(id));
            out.write(reinterpret_cast<const char*>(&max_layer_val), sizeof(max_layer_val));
            out.write(reinterpret_cast<const char*>(&max_layer_i_val), sizeof(max_layer_i_val));

            size_t total_size = (max_layer_val + 1) * (max_layer_i_val + 2) * sizeLinkList;

            auto ptr = linklist[id];
            out.write(ptr, total_size);
        }

        out.close();
    }

    static TrieIndex* loadIndex(const string& filename) {
        ifstream in(filename, ios::binary);
        if (!in) {
            throw runtime_error("Cannot open file for reading: " + filename);
        }

        // 1. Read basic parameters
        int M, ef_construction, dim, label_num;
        size_t maxNum, eleCount;
        double mult_, revSize_;

        in.read(reinterpret_cast<char*>(&M), sizeof(M));
        in.read(reinterpret_cast<char*>(&ef_construction), sizeof(ef_construction));
        in.read(reinterpret_cast<char*>(&dim), sizeof(dim));
        in.read(reinterpret_cast<char*>(&maxNum), sizeof(maxNum));
        in.read(reinterpret_cast<char*>(&eleCount), sizeof(eleCount));
        in.read(reinterpret_cast<char*>(&mult_), sizeof(mult_));
        in.read(reinterpret_cast<char*>(&revSize_), sizeof(revSize_));
        in.read(reinterpret_cast<char*>(&label_num), sizeof(label_num));

        int random_seed = 100;

        // 2. Create empty index object
        TrieIndex* index = new TrieIndex();
        index->level_generator_.seed(random_seed);
        index->M = M;
        index->ef_construction = ef_construction;
        index->dim = dim;
        index->maxNum = maxNum;
        index->eleCount = eleCount;
        index->mult_ = mult_;
        index->revSize_ = revSize_;
        index->label_num = label_num;
        index->max_label = label_num;

        // 3. Initialize space and distance functions
        index->space = hnswlib::L2Space(dim);
        index->fstdistfunc_ = index->space.get_dist_func();
        index->dist_func_param_ = index->space.get_dist_func_param();
        index->data_size_ = index->space.get_data_size();
        index->sizeLinkList = (M * sizeof(tableint) + sizeof(linklistsizeint));


        auto preMem = get_memory_usage();

        // 4. Read raw data arrays
        index->vecData_ = new char[maxNum * dim * sizeof(float)];
        in.read(index->vecData_, eleCount * dim * sizeof(float));

        auto postMem = get_memory_usage();
        print_memory_diff("vec", preMem, postMem);

        index->keyList_ = new int[maxNum];
        in.read(reinterpret_cast<char*>(index->keyList_), eleCount * sizeof(int));


        postMem = get_memory_usage();
        print_memory_diff("key", preMem, postMem);

        index->isDeleted = new bool[maxNum];
        in.read(reinterpret_cast<char*>(index->isDeleted), maxNum * sizeof(bool));


        postMem = get_memory_usage();
        print_memory_diff("del", preMem, postMem);

        // 5. Read item labels
        size_t num_items;
        in.read(reinterpret_cast<char*>(&num_items), sizeof(num_items));
        index->item_labels.resize(maxNum);
        for (size_t i = 0; i < num_items; i++) {
            size_t size;
            in.read(reinterpret_cast<char*>(&size), sizeof(size));
            index->item_labels[i].resize(size);
            in.read(reinterpret_cast<char*>(index->item_labels[i].data()),
                    size * sizeof(int));
        }


        postMem = get_memory_usage();
        print_memory_diff("label", preMem, postMem);

        index->eng_label = std::mt19937(37);
        index->hashedlabel.resize(label_num + 1);
        index->label_node.resize(label_num + 1);
        index->labelhash_generator_ = std::uniform_int_distribution<int>(1, label_num);
        for (int i = 1; i <= label_num; i++) {
            if (i <= index->label_the) index->hashedlabel[i] = i-1;
            else index->hashedlabel[i] = index->labelhash_generator_(index->eng_label) % (64-index->label_the) + index->label_the;
        }

        // max_layer, max_layer_i, labelhash
        index->max_layer.resize(maxNum);
        in.read(reinterpret_cast<char*>(index->max_layer.data()), maxNum * sizeof(int));

        index->max_layer_i.resize(maxNum);
        in.read(reinterpret_cast<char*>(index->max_layer_i.data()), maxNum * sizeof(int));

        index->labelhash.resize(maxNum);
        in.read(reinterpret_cast<char*>(index->labelhash.data()), maxNum * sizeof(unsigned long long));
        for(int i = 0; i < eleCount; i++) index->labelhash[i] = index->setBits(index->item_labels[i]);


        postMem = get_memory_usage();
        print_memory_diff("layer and hash", preMem, postMem);

        // 6. Initialize visited array
        index->visited_array = new unsigned int[maxNum]();
        index->tag = 0;
        index->tag2 = 0;

        // 7. Load Trie tree structure
        index->root = loadTrieNode(in, nullptr, index);


        postMem = get_memory_usage();
        print_memory_diff("trie", preMem, postMem);

        // 8. Load key2Id mapping
        size_t map_size;
        in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
        for (size_t i = 0; i < map_size; i++) {
            int key, id;
            in.read(reinterpret_cast<char*>(&key), sizeof(key));
            in.read(reinterpret_cast<char*>(&id), sizeof(id));
            index->key2Id[key] = id;
        }

        // 9. Load link lists
        index->linklist.resize(maxNum, nullptr);
        index->belong.resize(maxNum);

        size_t size_sum = 0;
        size_t ilayer_sum = 0;
        size_t layer_sum = 0;

        for (size_t i = 0; i < eleCount; i++) {
            size_t id;
            in.read(reinterpret_cast<char*>(&id), sizeof(id));
            int max_layer_val, max_layer_i_val;
            in.read(reinterpret_cast<char*>(&max_layer_val), sizeof(max_layer_val));
            index->max_layer[id] = max_layer_val;
            in.read(reinterpret_cast<char*>(&max_layer_i_val), sizeof(max_layer_i_val));
            index->max_layer_i[id] = max_layer_i_val;

            size_t total_size = (max_layer_val + 1) * (max_layer_i_val + 2) * index->sizeLinkList;
            size_sum += total_size;
            ilayer_sum += (max_layer_val + 1) * (max_layer_i_val + 2);
            layer_sum += max_layer_val + 1;
            char* ptr = static_cast<char*>(malloc(total_size));
            index->linklist[id] = ptr;
            in.read(ptr, total_size);
        }

        cout<<size_sum<<' '<<layer_sum<<' '<<ilayer_sum<<endl;
        postMem = get_memory_usage();
        print_memory_diff("graph", preMem, postMem);

        for(int i = 0 ;i < eleCount; i++){
            index->getbelong(i,index->item_labels[i]);
        }
        in.close();

        // 10. Preprocess label sets
        index->trie_label_load(index->root);
        std::cout << "Index loaded from " << filename << std::endl;
        return index;
    }


private:
public:
    TrieNode *root = nullptr;
    vector<bool> high_node;
    vector< vector<int > > item_labels;
    vector<int> max_layer;
    vector<int> max_layer_i;
    map<int,int> label_map;
    vector<pair<int,int> > sorted_array;
    vector<LabelType > label_list;
    vector<vector<TrieNode *> > label_node;

    vector<unsigned long long> labelhash;
    vector<int> hashedlabel;

    hnswlib::L2Space space;
    DISTFUNC<float> fstdistfunc_;
    void *dist_func_param_{nullptr};
    size_t data_size_{0};
    std::mt19937 eng;
    std::mt19937 eng_label;
    std::default_random_engine level_generator_;
    std::uniform_int_distribution<int> labelhash_generator_;
    double mult_{0.0}, revSize_{0.0};

    char* vecData_;
    int* keyList_;
    int* valueList_;
    bool* isDeleted;
    size_t maxNum, eleCount;
    size_t sizeLinkList;
    int M,M0;
    int skipLayer = 1;
    int ef_construction;
    int dim;
    double alpha = 1;
    int label_num;
    int label_the = 16;
    int max_label = 0;

    std::vector<char *> linklist;

    unsigned int *visited_array;
    unsigned int tag = 0, tag2 = 0;
    std::vector<int> sortedArray;
    vector<vector<TrieNode *> >belong;

    std::unordered_map<int,int> key2Id;

    int numEdges = 0;

    unsigned long long setBits(const std::vector<int>& a) {
        unsigned long long result = 0;
        if (label_num <= label_the) {
            for (int pos : a) {
                if (pos > label_num) continue;
                if (pos < 1 || pos > 64) {
                    throw std::out_of_range("Position must be between 1 and 64");
                }
                result |= (1ULL << (pos - 1));
            }
        }
        else {
            for (int pos : a) {
                if (pos > label_num) continue;
                result |= (1ULL << hashedlabel[pos]);
            }
        }
        return result;
    }

    int getRandomLevel(double reverse_size) {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double r = -log(distribution(level_generator_)) * reverse_size;
        return (int) r;
    }

    double countR = 0;

    void build_index(TrieNode *node, size_t begin, size_t end, int idx, int layer){
        node->size = node->actual_size = node->graph_size = end - begin + 1;
        node->other_size = 0;
        if(alpha == 1 || layer== 0 || ( (int)(log(node->size)/log(alpha)) < (int)(log(node->parent->pre_index->size) / log(alpha)))){
            tableint ep_id = -1;
            node->layer = layer;
            node->pre_index = node;
            for(int i = begin; i <= end; i++) {

                int id = sorted_array[i].second;

                belong[id].resize(layer + 1);
                belong[id][layer] = node;
                max_layer[id] = layer;
                tableint epp = ep_id;
                for(int j = max_layer_i[id]; j >=0 ; j--) {
                    if(epp == -1 || max_layer_i[epp] < j) {
                        linklistsizeint *ll_cur = get_linklist(id, layer,j);

                        setListCount(ll_cur, 0);
                        continue;
                    }
                    std::vector<tableint> ep_ids = {epp};
                    auto candidates = searchBaseLayer(ep_ids, getDataByInternalId(id), layer,j);

                    if(candidates.size() > 0) {
                        if(j ==0 )getNeighborsByHeuristic2(candidates, 2 * M);
                        else getNeighborsByHeuristic2(candidates, M);
                        epp = connectEdges(getDataByInternalId(id), id, candidates, layer, j);
                        int *data = (int *) get_linklist(id, layer, j);
                        size_t size = getListCount((linklistsizeint *) data);
                    }
                }
                if(ep_id == -1 || max_layer_i[ep_id] < max_layer_i[id]){
                    ep_id = id;
                }
            }
            node->entry = ep_id;
            layer ++;
        }
        else{
            node->pre_index = node->parent->pre_index;
            node->parent->son = node;
        }
        if(item_labels[sorted_array[begin].second].size() == idx){
            return;
        }

        while(item_labels[sorted_array[end].second].size() == idx){
            end--;
        }

        for(int i = begin; i <= end; i++){
            sorted_array[i].first = item_labels[sorted_array[i].second][idx];
        }

        sort(sorted_array.begin()+(int)begin,sorted_array.begin()+(int)end + 1);



        size_t pre = begin;
        for(int i = begin + 1; i<= end; i++){
            if(item_labels[sorted_array[i].second][idx]!=item_labels[sorted_array[pre].second][idx]){
                build_index(node->child[item_labels[sorted_array[pre].second][idx]],pre,i-1,idx + 1, layer);
                pre = i;
            }
        }
        build_index(node->child[item_labels[sorted_array[pre].second][idx]],pre,end,idx + 1, layer);
    }

    void adjustTrie(int id, const vector<int> &label_set){
        TrieNode *node = root;
        if(root == nullptr) node = root = new TrieNode(nullptr, 0);
        for (int i = 0; i < label_set.size(); i++) {
            if(node->pre_index == node){
                if(int(log(node->size)/log(alpha)) > int(log(node->size-1)/log(alpha))){

                    vector<int> node_id_set;
                    collectELments(node, node_id_set);

                    if(node->son == nullptr){
                        if((node != root) && (int(log(node->size)/log(alpha)) == int(log(node->parent->pre_index->size)/log(alpha)))){
                            node->pre_index = node->parent->pre_index;
                            node->parent->son = node;
                            for(auto eid: node_id_set){
                                belong[eid][node->layer] = nullptr;
                            }
                        }
                    }
                    else{

                        vector<int> son_id_set;
                        collectELments(node->son, son_id_set);

                        if((node != root) && (int(log(node->size)/log(alpha)) == int(log(node->parent->pre_index->size)/log(alpha)))){
                            node->pre_index = node->parent->pre_index;
                            node->parent->son = node;
                            for(auto eid: node_id_set){
                                belong[eid][node->layer] = nullptr;
                            }

                            node->son->layer = node->layer;
                            node->son->pre_index = node->son;
                            int ep_id = -1;
                            for(auto eid: son_id_set){
                                if(eid == id) continue;
                                if(ep_id == -1 || max_layer_i[ep_id] < max_layer_i[eid]){
                                    ep_id = eid;
                                }
                            }
                            node->son->entry = ep_id;
                            node->son->graph_size = node->graph_size;
                            auto tmpnode = node->son;
                            while (tmpnode->son!=nullptr)
                            {
                                tmpnode = tmpnode->son;
                                tmpnode->pre_index = node->son;
                            }
                            for(auto eid: son_id_set){
                                belong[eid][node->layer] = node->son;
                            }
                        }
                        else{
                            for(auto eid: son_id_set){
                                max_layer[eid] += 1;
                                if(eid == id) continue;
                                linklist[eid] = (char *) realloc(linklist[eid], (max_layer[eid] + 1) * sizeLinkList * (max_layer_i[eid] + 2));
                            }

                            changeLayer(node->son,1);
                            for(auto eid: son_id_set){
                                belong[eid].resize(max_layer[eid] + 1);
                                for(int l = max_layer[eid]; l > node->layer; l--){
                                    belong[eid][l] = belong[eid][l-1];
                                    if(eid == id) continue;
                                    auto ptr_old = get_linklist(eid, l-1, 0);
                                    auto ptr_new = get_linklist(eid, l, 0);
                                    memcpy(ptr_new, ptr_old, sizeLinkList * (max_layer_i[eid] + 2));
                                }
                            }
                            node->son->layer = node->layer + 1;
                            node->son->pre_index = node->son;
                            int ep_id = -1;
                            for(auto eid: son_id_set){
                                if(eid == id) continue;
                                if(ep_id == -1 || max_layer_i[ep_id] < max_layer_i[eid]){
                                    ep_id = eid;
                                }
                            }
                            node->son->entry = ep_id;

                            node->son->graph_size = node->graph_size;
                            auto tmpnode = node->son;
                            while (tmpnode->son!=nullptr)
                            {
                                tmpnode = tmpnode->son;
                                tmpnode->pre_index = node->son;
                            }
                            for(auto eid: son_id_set){
                                belong[eid][node->son->layer] = node->son;
                            }
                        }

                        if(log(node->graph_size)/log(alpha) > int(log(node->son->size)/log(alpha))+1){
                            rebuild(node->son, son_id_set,id);
                        }

                        node->son = nullptr;
                    }
                }
            }
            node = node->child[label_set[i]];
        }
    }

    void rebuild(TrieNode *node, const vector<int> &id_set, int no_id = -1){
        node->graph_size = id_set.size();
        int layer = node->layer;
        tableint ep_id = -1;
        node->pre_index = node;
        for(int i = 0; i < id_set.size(); i++) {

            int id = id_set[i];
            if(id == no_id) continue;

            if(belong[id].size()<layer+1) belong[id].resize(layer + 1);
            belong[id][layer] = node;
            tableint epp = ep_id;
            for(int j = max_layer_i[id]; j >=0 ; j--) {
                if(epp == -1 || max_layer_i[epp] < j) {
                    linklistsizeint *ll_cur = get_linklist(id, layer,j);

                    setListCount(ll_cur, 0);
                    continue;
                }
                std::vector<tableint> ep_ids = {epp};
                auto candidates = searchBaseLayer(ep_ids, getDataByInternalId(id), layer,j);

                if(candidates.size() > 0) {
                    if(j ==0 )getNeighborsByHeuristic2(candidates, 2 * M);
                    else getNeighborsByHeuristic2(candidates, M);
                    epp = connectEdges2(getDataByInternalId(id), id, candidates, layer, j);
                }
            }
            if(ep_id == -1 || max_layer_i[ep_id] < max_layer_i[id]){
                ep_id = id;
            }
        }
        node->entry = ep_id;

    }

    void changeLayer(TrieNode *node, int delta){
        if(node->child.empty()) return;
        for(auto ptr: node->child){
            ptr.second->layer += delta;
            changeLayer(ptr.second, delta);
        }
    }

    void deleteFromTrie(int id, const vector<int> &label_set){
        TrieNode *node = root;
        for (int i = 0; i < label_set.size(); i++) {
            node->actual_size--;
            if(node->actual_size < node->graph_size/2){
                node->graph_size = node->actual_size;
                vector<int> id_set;
                collectELments(node, id_set);
                rebuild(node, id_set,id);
            }
            node = node->child[label_set[i]];
        }
        for(int idx = 0; idx < node->ids.size(); idx++){
            if(node->ids[idx] == id){
                node->ids.erase(node->ids.begin() + idx);
                break;
            }
        }
    }

    void collectELments(TrieNode *node, vector<int> &id_set){
        if(node->child.empty()){
            for(auto id: node->ids){
                id_set.push_back(id);
            }
            return;
        }
        for(auto ptr: node->child){
            collectELments(ptr.second, id_set);
        }
    }

    tableint insertIntoGraph(TrieNode *node, int id, const vector<int> &label_set, char* data,int idx){
        tableint bottom_ep_id = -1;
        if(node->child.empty()){
        }
        else{
            bottom_ep_id = insertIntoGraph(node->child[label_set[idx]], id, label_set, data, idx + 1);
        }
        if(node->pre_index == node){
            belong[id][node->layer] = node;
            if(node->entry == -1){
                node->entry = id;
                for(int j = max_layer_i[id]; j >=0 ; j--) {
                    linklistsizeint *ll_cur = get_linklist(id, node->layer,j);
                    setListCount(ll_cur, 0);
                }
                return -1;
            }

            tableint ep_id = node->entry;
            for(int j = max_layer_i[id]; j >=0 ; j--) {
                if(ep_id == -1||max_layer_i[node->entry] < j) {
                    linklistsizeint *ll_cur = get_linklist(id, node->layer,j);

                    setListCount(ll_cur, 0);
                    continue;
                }
                if(j == 0 && bottom_ep_id != -1){
                    ep_id = bottom_ep_id;
                }
                std::vector<tableint> ep_ids = {ep_id};
                auto candidates = searchBaseLayer(ep_ids, getDataByInternalId(id), node->layer,j);

                if(candidates.size() > 0) {
                    if(j ==0 )getNeighborsByHeuristic2(candidates, 2 * M);
                    else getNeighborsByHeuristic2(candidates, M);
                    ep_id = connectEdges2(getDataByInternalId(id), id, candidates, node->layer, j);
                    bottom_ep_id = ep_id;
                }
            }
            if(node->entry == -1|| max_layer_i[node->entry] < max_layer_i[id]){
                node->entry = id;
            }
        }
        return bottom_ep_id;
    }


    bool superset_search(TrieNode *root,
                     const std::vector<int> &label_set,
                     int /*idx*/,
                     std::vector<TrieNode *> &node_set) {
    if (label_set.empty()) return false;

    std::unordered_map<TrieNode*, TrieNode*> first_hit_child;
    std::unordered_set<TrieNode*> pushed;

    bool any_hit = false;
    int last = label_set.back();

    if (last < 0 || (size_t)last >= label_node.size()) {
        return false;
    }

    for (TrieNode* leaf : label_node[last]) {
        TrieNode* cur = leaf;
        int pos = (int)label_set.size() - 1;

        while (cur != nullptr && pos >= 0) {
            if (cur->parent == nullptr) break;
            auto it = cur->parent->child.find(label_set[pos]);
            if (it != cur->parent->child.end() && it->second == cur) {
                pos--;
            }
            cur = cur->parent;
        }

        if (pos >= 0) continue;

        any_hit = true;

        TrieNode* out = leaf->pre_index;
        if (pushed.insert(out).second) {
            node_set.push_back(out);
        }
        leaf->pre_index->flag = tag2;

        TrieNode* child_on_path = leaf;
        TrieNode* u = leaf->parent;
        while (u != nullptr) {
            auto it = first_hit_child.find(u);
            if (it == first_hit_child.end()) {
                first_hit_child.emplace(u, child_on_path);
            } else if (it->second != child_on_path) {
                u->pre_index->flag = tag2;
            }
            child_on_path = u;
            u = u->parent;
        }
    }

    return any_hit;
}



    int superset_search(TrieNode *node, const vector<int> &label_set, int idx) {
        if (idx == label_set.size()) {
            return 1;
        }
        bool flag = false;
        int from = idx > 0 ? label_set[idx - 1] : 0;
        int upto = label_set[idx];
        int cnt = 0;
        for (auto ptr: node->child) {
            bool exist = false;
            if(ptr.first > from && ptr.first < upto){
                 cnt += superset_search(ptr.second,label_set,idx);
            } else if (ptr.first == upto){
                cnt += superset_search(ptr.second,label_set,idx + 1);
            }
        }
        return cnt;
    }

    bool subset_search(TrieNode *node, const vector<int> &label_set, int idx,
                         vector<TrieNode *> &node_set) {
        if (node->child.empty()) {
            node_set.push_back(node);
            node->flag = true;
            return true;
        }
        if (idx == label_set.size()) {
            return false;
        }
        bool flag = false;
        for (auto ptr: node->child) {
            bool exist = false;
            if(ptr.first == label_set[idx]){
                exist = subset_search(ptr.second,label_set,idx + 1,node_set);
            } else if (ptr.first > label_set[idx] && idx < label_set.size()){
                exist = subset_search(node,label_set,idx + 1,node_set);
            }
            flag = exist || flag;
            if(exist && flag) node->flag = true;
        }

        return flag;
    }

    void equal_search(TrieNode *node, const vector<int> &label_set, int idx,
                       vector<TrieNode *> &node_set) {
        if (idx == label_set.size() && node->child.count(label_num + 1) != 0) {
            node_set.push_back(node->child[label_num + 1]->pre_index);
            node->child[label_num + 1]->pre_index->flag = tag2;
            return;
        }
        if(node->child.count(label_set[idx]) != 0)
            equal_search(node->child[label_set[idx]],label_set,idx+1,node_set);
        else {
            std::cout << "No equal labels found." << std::endl;
        }
    }

    bool overlap_search(
    TrieNode *root,
    const std::vector<int> &label_set,
    int idx,
    std::vector<TrieNode *> &node_set)
{
    if (label_set.empty()) return false;

    std::unordered_map<TrieNode*, TrieNode*> first_hit_child;
    std::unordered_set<TrieNode*> pushed;
    std::unordered_set<TrieNode*> seen_leaf;

    bool any_hit = false;

    for (int lab : label_set) {
        if (lab < 0 || (size_t)lab >= label_node.size()) continue;

        auto &vec = label_node[lab];
        if (vec.empty()) continue;

        for (TrieNode* leaf : vec) {
            if (!leaf) continue;
            if (!seen_leaf.insert(leaf).second) continue;

            any_hit = true;

            TrieNode* out = leaf->pre_index;
            if (out && pushed.insert(out).second) {
                node_set.push_back(out);
            }

            if (leaf->pre_index) leaf->pre_index->flag = tag2;

            TrieNode* child_on_path = leaf;
            TrieNode* u = leaf->parent;
            while (u != nullptr) {
                auto it = first_hit_child.find(u);
                if (it == first_hit_child.end()) {
                    first_hit_child.emplace(u, child_on_path);
                } else if (it->second != child_on_path) {
                    if (u->pre_index) u->pre_index->flag = tag2;
                }
                child_on_path = u;
                u = u->parent;
            }
        }
    }

    return any_hit;
}


    ResultHeap searchBaseLayer(const std::vector<tableint> &ep_ids, const void *data_point, int layer, int ilayer) {
        tag ++;

        ResultHeap top_candidates;
        ResultHeap candidateSet;

        float lowerBound;

        TrieNode *nd = belong[ep_ids[0]][layer];

        if(ep_ids[0] == -1) return {};

        for(int i = 0; i < ep_ids.size(); i++) {
            int ep_id = ep_ids[i];
            float dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
            if(!isDeleted[ep_id]) {
                top_candidates.emplace(dist, ep_id);
                candidateSet.emplace(-dist, ep_id);
            }
            else{
                candidateSet.emplace(-std::numeric_limits<float>::max(), ep_id);
            }
            visited_array[ep_id] = tag;
        }

        if(!top_candidates.empty())
            lowerBound = top_candidates.top().first;
        else
            lowerBound = std::numeric_limits<float>::max();


        while (!candidateSet.empty()) {
            std::pair<float, tableint> curr_el_pair = candidateSet.top();
            if ((-curr_el_pair.first) > lowerBound && top_candidates.size() == ef_construction) {
                break;
            }
            candidateSet.pop();

            tableint curNodeNum = curr_el_pair.second;

            int *data = (int *) get_linklist(curNodeNum, layer, ilayer);
            size_t size = getListCount((linklistsizeint *) data);
            tableint *datal = (tableint *) (data + 1);

            for (size_t j = 0; j < size; j++) {
                tableint candidate_id = *(datal + j);
                if(max_layer[candidate_id] < layer) continue;
                if(belong [candidate_id][layer] != nd) continue;
#ifdef USE_SSE
                _mm_prefetch((char *) (visited_array + *(datal + j + 1)), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
                _mm_prefetch((char *) (visited_array + *(datal + j + 2)), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*(datal + j + 2)), _MM_HINT_T0);
#endif
                if (visited_array[candidate_id] == tag) continue;
                visited_array[candidate_id] = tag;
                char *currObj1 = (getDataByInternalId(candidate_id));

                float dist1 = fstdistfunc_(data_point, currObj1, dist_func_param_);
                if (top_candidates.size() < ef_construction || lowerBound > dist1) {
                    candidateSet.emplace(-dist1, candidate_id);
#ifdef USE_SSE
                    _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif

                    if(!isDeleted[candidate_id])
                        top_candidates.emplace(dist1, candidate_id);
                    if (top_candidates.size() > ef_construction)
                        top_candidates.pop();

                    if (!top_candidates.empty())
                        lowerBound = top_candidates.top().first;
                }
            }
        }

        return top_candidates;
    }

    tableint
    findEntry(const void *query_data, TrieNode *nd, int lastLayer) const {
        tableint currObj = nd->entry;
        float curdist = fstdistfunc_(query_data, getDataByInternalId(currObj), dist_func_param_);

        int layer = nd->layer;

        for (int ilayer = max_layer_i[nd->entry]; ilayer > 0 ; ilayer --) {
            bool changed = true;
            while (changed) {
                changed = false;
                unsigned int *data;

                data = (unsigned int *) get_linklist(currObj, layer, ilayer);
                int size = getListCount(data);

                tableint *datal = (tableint *) (data + 1);
                for (int i = 0; i < size; i++) {
                    tableint cand = datal[i];
#ifdef USE_SSE
                    _mm_prefetch(getDataByInternalId(*(datal + i + 1)), _MM_HINT_T0);
                    _mm_prefetch(getDataByInternalId(*(datal + i + 2)), _MM_HINT_T0);
#endif
                    float d = fstdistfunc_(query_data, getDataByInternalId(cand), dist_func_param_);

                    if (d < curdist) {
                        curdist = d;
                        currObj = cand;
                        changed = true;
                    }
                }
            }
        }
        return currObj;
    }

    tableint connectEdges(
            const void *data_point,
            tableint cur_c,
            ResultHeap &top_candidates,
            int layer,
            int ilayer) {

        std::vector<tableint> selectedNeighbors;
        int M_ = M;
        if(ilayer ==0 ) M_ = 2 * M;
        selectedNeighbors.reserve(M_);
        while (top_candidates.size() > 0) {
            selectedNeighbors.push_back(top_candidates.top().second);
            top_candidates.pop();
        }

        tableint next_closest_entry_point = selectedNeighbors.back();

        {
            linklistsizeint *ll_cur = get_linklist(cur_c, layer,ilayer);

            setListCount(ll_cur, selectedNeighbors.size());
            tableint *data = (tableint *) (ll_cur + 1);
            for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {
                data[idx] = selectedNeighbors[idx];
            }
        }

        for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {

            linklistsizeint *ll_other = get_linklist(selectedNeighbors[idx], layer,ilayer);

            size_t sz_link_list_other = getListCount(ll_other);

            tableint *data = (tableint *) (ll_other + 1);
            if (sz_link_list_other < M_) {
                data[sz_link_list_other] = cur_c;
                setListCount(ll_other, sz_link_list_other + 1);
            } else {
                // finding the "weakest" element to replace it with the new one
                float d_max = fstdistfunc_(getDataByInternalId(cur_c), getDataByInternalId(selectedNeighbors[idx]),
                                           dist_func_param_);
                ResultHeap candidates;
                candidates.emplace(d_max, cur_c);

                for (size_t j = 0; j < sz_link_list_other; j++) {
                    candidates.emplace(
                            fstdistfunc_(getDataByInternalId(data[j]), getDataByInternalId(selectedNeighbors[idx]),
                                         dist_func_param_), data[j]);
                }

                if(ilayer ==0 )getNeighborsByHeuristic2(candidates, 2 * M);
                else getNeighborsByHeuristic2(candidates, M_);

                int indx = 0;
                while (candidates.size() > 0) {
                    data[indx] = candidates.top().second;
                    candidates.pop();
                    indx++;
                }

                setListCount(ll_other, indx);
            }
        }

        return next_closest_entry_point;
    }

    tableint connectEdges2(
            const void *data_point,
            tableint cur_c,
            ResultHeap &top_candidates,
            int layer,
            int ilayer) {

        std::vector<tableint> selectedNeighbors;
        selectedNeighbors.reserve(M);
        while (top_candidates.size() > 0) {
            selectedNeighbors.push_back(top_candidates.top().second);
            top_candidates.pop();
        }

        tableint next_closest_entry_point = selectedNeighbors.back();

        {
            linklistsizeint *ll_cur = get_linklist(cur_c, layer,ilayer);

            setListCount(ll_cur, selectedNeighbors.size());
            tableint *data = (tableint *) (ll_cur + 1);
            for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {
                data[idx] = selectedNeighbors[idx];
            }
        }

        for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {
            if(belong[selectedNeighbors[idx]][layer]!=belong[cur_c][layer]) continue;

            linklistsizeint *ll_other = get_linklist(selectedNeighbors[idx], layer,ilayer);

            size_t sz_link_list_other = getListCount(ll_other);

            tableint *data = (tableint *) (ll_other + 1);
            if (sz_link_list_other < M) {
                data[sz_link_list_other] = cur_c;
                setListCount(ll_other, sz_link_list_other + 1);
            } else {
                // finding the "weakest" element to replace it with the new one
                float d_max = fstdistfunc_(getDataByInternalId(cur_c), getDataByInternalId(selectedNeighbors[idx]),
                                           dist_func_param_);
                ResultHeap candidates;
                candidates.emplace(d_max, cur_c);

                for (size_t j = 0; j < sz_link_list_other; j++) {
                    candidates.emplace(
                            fstdistfunc_(getDataByInternalId(data[j]), getDataByInternalId(selectedNeighbors[idx]),
                                         dist_func_param_), data[j]);
                }

                if(ilayer ==0 )getNeighborsByHeuristic2(candidates, 2 * M);
                else getNeighborsByHeuristic2(candidates, M);

                int indx = 0;
                while (candidates.size() > 0) {
                    data[indx] = candidates.top().second;
                    candidates.pop();
                    indx++;
                }

                setListCount(ll_other, indx);
            }
        }

        return next_closest_entry_point;
    }

    void getNeighborsByHeuristic2(
            ResultHeap &top_candidates,
            const size_t M) {
        if (top_candidates.size() < M) {
            return;
        }

        std::priority_queue<std::pair<float, tableint>> queue_closest;
        std::vector<std::pair<float, tableint>> return_list;
        while (top_candidates.size() > 0) {
            queue_closest.emplace(-top_candidates.top().first, top_candidates.top().second);
            top_candidates.pop();
        }

        while (queue_closest.size()) {
            if (return_list.size() >= M)
                break;
            std::pair<float, tableint> curent_pair = queue_closest.top();
            float dist_to_query = -curent_pair.first;
            queue_closest.pop();
            bool good = true;

            for (std::pair<float, tableint> second_pair : return_list) {
                float curdist =
                        fstdistfunc_(getDataByInternalId(second_pair.second),
                                     getDataByInternalId(curent_pair.second),
                                     dist_func_param_);
                if (curdist < dist_to_query) {
                    good = false;
                    break;
                }
            }
            if (good) {
                return_list.push_back(curent_pair);
            }
        }

        for (std::pair<float, tableint> curent_pair : return_list) {
            top_candidates.emplace(-curent_pair.first, curent_pair.second);
        }
    }

    int comp;
    double time;

    inline bool isContained(unsigned long long qhash, vector<int> &ql, int id) {
        if ((qhash & labelhash[id]) != qhash) return false;
        if (!isContainment(ql, item_labels[id])) return false;
        return true;
    }

    inline bool isOverlaped(unsigned long long qhash, vector<int> &ql, int id) {
        if ((qhash & labelhash[id]) == 0) return false;
        if (!isOverlap(ql, item_labels[id])) return false;
        return true;
    }

    ResultHeap
    searchBaseLayer0(std::vector<tableint> ep_ids, const void *data_point, unsigned long long query_hash, vector<int> &ql, int ef,int layer=0) {

        ResultHeap top_candidates;
        ResultHeap candidateSet;

        float lowerBound;
        for(int i = 0; i < ep_ids.size(); i++) {
            int ep_id = ep_ids[i];
            if(visited_array[ep_id] == tag) continue;
            float dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
            if(!isDeleted[ep_id] && isContained(query_hash, ql, ep_id)) {
                top_candidates.emplace(dist, ep_id);
                candidateSet.emplace(-dist, ep_id);
            }
            else{
                candidateSet.emplace(-std::numeric_limits<float>::max(), ep_id);
            }
            visited_array[ep_id] = tag;
        }
        if(!top_candidates.empty())
            lowerBound = top_candidates.top().first;
        else
            lowerBound = std::numeric_limits<float>::max();
        while (!candidateSet.empty()) {
            std::pair<float, tableint> curr_el_pair = candidateSet.top();
            tableint curNodeNum = curr_el_pair.second;
            if ((-curr_el_pair.first) > lowerBound && top_candidates.size() == ef) {
                break;
            }
            candidateSet.pop();

            int cntM = 0;
            for(int i = 0; i <= max_layer[curNodeNum]; i++) {
                if(belong[curNodeNum][i] == nullptr) continue;
                if(!(belong[curNodeNum][i]->flag==tag2)) continue;
                int *data = (int *) get_linklist(curNodeNum, i, 0);


                size_t size = getListCount((linklistsizeint *) data);
                tableint *datal = (tableint *) (data + 1);
#ifdef USE_SSE
                _mm_prefetch((char *) (visited_array + *(data + 1)), _MM_HINT_T0);
                _mm_prefetch((char *) (visited_array + *(data + 1) + 64), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*datal), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*(datal + 1)), _MM_HINT_T0);
#endif

                for (size_t j = 0; j < size; j++) {
                    tableint candidate_id = *(datal + j);
#ifdef USE_SSE
                        _mm_prefetch((char *) (visited_array + *(datal + j + 1)), _MM_HINT_T0);
                        _mm_prefetch(getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
                        _mm_prefetch((char *) (visited_array + *(datal + j + 2)), _MM_HINT_T0);
                        _mm_prefetch(getDataByInternalId(*(datal + j + 2)), _MM_HINT_T0);
#endif
                    if (visited_array[candidate_id] == tag) continue;
                    visited_array[candidate_id] = tag;
                    char *currObj1 = (getDataByInternalId(candidate_id));

                    tableint cid = candidate_id;

                    if(!isContained(query_hash, ql, cid)) {
                        continue;
                    }
                    cntM ++;
                    float dist1 = fstdistfunc_(data_point, currObj1, dist_func_param_);
                    if (top_candidates.size() < ef || lowerBound > dist1) {
                        candidateSet.emplace(-dist1, cid);
#ifdef USE_SSE
                        _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif


                        if(!isDeleted[candidate_id])
                            top_candidates.emplace(dist1, cid);
                        if (top_candidates.size() > ef)
                            top_candidates.pop();

                        if (!top_candidates.empty())
                            lowerBound = top_candidates.top().first;
                    }
                }
                if(cntM > 2*M)
                    break;
            }

        }

        return top_candidates;
    }

    ResultHeap
    searchBaseLayer0forEquality(std::vector<tableint> ep_ids, const void *data_point, unsigned long long query_hash, vector<int> &ql, int ef,int layer=0) {

        ResultHeap top_candidates;
        ResultHeap candidateSet;

        float lowerBound;
        for(int i = 0; i < ep_ids.size(); i++) {
            int ep_id = ep_ids[i];
            if(visited_array[ep_id] == tag) continue;
            float dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
            if(!isDeleted[ep_id]) {
                top_candidates.emplace(dist, ep_id);
                candidateSet.emplace(-dist, ep_id);
            }
            else{
                candidateSet.emplace(-std::numeric_limits<float>::max(), ep_id);
            }
            visited_array[ep_id] = tag;
        }
        if(!top_candidates.empty())
            lowerBound = top_candidates.top().first;
        else
            lowerBound = std::numeric_limits<float>::max();
        while (!candidateSet.empty()) {
            std::pair<float, tableint> curr_el_pair = candidateSet.top();
            tableint curNodeNum = curr_el_pair.second;
            if ((-curr_el_pair.first) > lowerBound && top_candidates.size() == ef) {
                break;
            }
            candidateSet.pop();

            int cntM = 0;
            for(int i = 0; i <= max_layer[curNodeNum]; i++) {
                if(belong[curNodeNum][i] == nullptr) continue;
                if(!(belong[curNodeNum][i]->flag==tag2)) continue;
                int *data = (int *) get_linklist(curNodeNum, i, 0);
                size_t size = getListCount((linklistsizeint *) data);
                tableint *datal = (tableint *) (data + 1);
#ifdef USE_SSE
                _mm_prefetch((char *) (visited_array + *(data + 1)), _MM_HINT_T0);
                _mm_prefetch((char *) (visited_array + *(data + 1) + 64), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*datal), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*(datal + 1)), _MM_HINT_T0);
#endif

                for (size_t j = 0; j < size; j++) {
                    tableint candidate_id = *(datal + j);
#ifdef USE_SSE
                        _mm_prefetch((char *) (visited_array + *(datal + j + 1)), _MM_HINT_T0);
                        _mm_prefetch(getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
                        _mm_prefetch((char *) (visited_array + *(datal + j + 2)), _MM_HINT_T0);
                        _mm_prefetch(getDataByInternalId(*(datal + j + 2)), _MM_HINT_T0);
#endif
                    if (visited_array[candidate_id] == tag) continue;
                    visited_array[candidate_id] = tag;
                    char *currObj1 = (getDataByInternalId(candidate_id));

                    tableint cid = candidate_id;
                    cntM ++;
                    float dist1 = fstdistfunc_(data_point, currObj1, dist_func_param_);


                    if (top_candidates.size() < ef || lowerBound > dist1) {
                        candidateSet.emplace(-dist1, cid);

#ifdef USE_SSE
                        _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif


                        if(!isDeleted[candidate_id])
                            top_candidates.emplace(dist1, cid);
                        if (top_candidates.size() > ef)
                            top_candidates.pop();

                        if (!top_candidates.empty())
                            lowerBound = top_candidates.top().first;
                    }
                }
                if(cntM > 2*M)
                    break;
            }

        }

        return top_candidates;
    }


    ResultHeap
    searchBaseLayer0ForOverlap(std::vector<tableint> ep_ids, const void *data_point, unsigned long long query_hash, vector<int> &ql, int ef) {

        ResultHeap top_candidates;
        ResultHeap candidateSet;

        float lowerBound;
        for(int i = 0; i < ep_ids.size(); i++) {
            int ep_id = ep_ids[i];
            if(visited_array[ep_id] == tag) continue;
            float dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
            if(!isDeleted[ep_id] && isOverlaped(query_hash, ql, ep_id)) {
                top_candidates.emplace(dist, ep_id);
                candidateSet.emplace(-dist, ep_id);
            }
            else{
                candidateSet.emplace(-std::numeric_limits<float>::max(), ep_id);
            }
            visited_array[ep_id] = tag;
        }

        if(!top_candidates.empty())
            lowerBound = top_candidates.top().first;
        else
            lowerBound = std::numeric_limits<float>::max();
        while (!candidateSet.empty()) {
            std::pair<float, tableint> curr_el_pair = candidateSet.top();
            tableint curNodeNum = curr_el_pair.second;
            if ((-curr_el_pair.first) > lowerBound && top_candidates.size() == ef) {
                break;
            }
            candidateSet.pop();

            int cntM = 0;
            for(int i = max_layer[curNodeNum]; i >= 0; i--) {
                if(belong[curNodeNum][i] == nullptr) continue;
                if(!(belong[curNodeNum][i]->flag==tag2)) continue;
                int *data = (int *) get_linklist(curNodeNum, i, 0);


                size_t size = getListCount((linklistsizeint *) data);
                tableint *datal = (tableint *) (data + 1);
#ifdef USE_SSE
                _mm_prefetch((char *) (visited_array + *(data + 1)), _MM_HINT_T0);
                _mm_prefetch((char *) (visited_array + *(data + 1) + 64), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*datal), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*(datal + 1)), _MM_HINT_T0);
#endif

                for (size_t j = 0; j < size; j++) {
                    tableint candidate_id = *(datal + j);
#ifdef USE_SSE
                        _mm_prefetch((char *) (visited_array + *(datal + j + 1)), _MM_HINT_T0);
                        _mm_prefetch(getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
                        _mm_prefetch((char *) (visited_array + *(datal + j + 2)), _MM_HINT_T0);
                        _mm_prefetch(getDataByInternalId(*(datal + j + 2)), _MM_HINT_T0);
#endif
                    if (visited_array[candidate_id] == tag) continue;
                    visited_array[candidate_id] = tag;
                    char *currObj1 = (getDataByInternalId(candidate_id));

                    tableint cid = candidate_id;

                    if(!isOverlaped(query_hash, ql, cid)) {
                        continue;
                    }
                    cntM ++;
                    float dist1 = fstdistfunc_(data_point, currObj1, dist_func_param_);
                    if (top_candidates.size() < ef || lowerBound > dist1) {
                            candidateSet.emplace(-dist1, cid);
#ifdef USE_SSE
                        _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif


                        if(!isDeleted[candidate_id])
                            top_candidates.emplace(dist1, cid);
                        if (top_candidates.size() > ef)
                            top_candidates.pop();

                        if (!top_candidates.empty())
                            lowerBound = top_candidates.top().first;
                    }
                }
                if(cntM > 2*M)
                    break;
            }

        }

        return top_candidates;
    }


    inline char *getDataByInternalId(tableint internal_id) const {
        return (char*)(vecData_ + internal_id * data_size_);
    }

    linklistsizeint *get_linklist(tableint internal_id, int layer, int iLayer) const {
        if (iLayer == 0)  return (linklistsizeint *) (linklist[internal_id] + sizeLinkList * layer * (max_layer_i[internal_id] + 2));
        return (linklistsizeint *) (linklist[internal_id] + sizeLinkList * layer * (max_layer_i[internal_id] + 2) + sizeLinkList * (iLayer + 1));
    }


    unsigned short int getListCount(linklistsizeint * ptr) const {
        return *((unsigned short int *)ptr);
    }

    void setListCount(linklistsizeint * ptr, unsigned short int size) const {
        *((unsigned short int*)(ptr))=*((unsigned short int *)&size);
    }

    void saveTrieNode(ofstream& out, TrieNode* node) {
        if (!node) {
            int null_marker = -1;
            out.write(reinterpret_cast<const char*>(&null_marker), sizeof(null_marker));
            return;
        }


        int null_marker = 1;
        out.write(reinterpret_cast<const char*>(&null_marker), sizeof(null_marker));
        bool has_index = node->pre_index == node;

        // Save node metadata
        out.write(reinterpret_cast<const char*>(&node->layer), sizeof(node->layer));
        out.write(reinterpret_cast<const char*>(&node->entry), sizeof(node->entry));
        out.write(reinterpret_cast<const char*>(&node->flag), sizeof(node->flag));
        out.write(reinterpret_cast<const char*>(&node->size), sizeof(node->size));
        out.write(reinterpret_cast<const char*>(&node->other_size), sizeof(node->other_size));
        out.write(reinterpret_cast<const char*>(&node->actual_size), sizeof(node->actual_size));
        out.write(reinterpret_cast<const char*>(&node->graph_size), sizeof(node->graph_size));
        out.write(reinterpret_cast<const char*>(&has_index), sizeof(has_index));

        // Save element IDs
        size_t element_count = node->ids.size();
        out.write(reinterpret_cast<const char*>(&element_count), sizeof(element_count));
        for(auto id:node->ids){
            out.write(reinterpret_cast<const char*>(&id), sizeof(id));
        }

        // Save children
        size_t child_count = node->child.size();
        out.write(reinterpret_cast<const char*>(&child_count), sizeof(child_count));
        for (const auto& child_pair : node->child) {
            int key = child_pair.first;
            out.write(reinterpret_cast<const char*>(&key), sizeof(key));
            saveTrieNode(out, child_pair.second);
        }
    }

    static int getlog(int a, int b){
        return (int)(log(a)/log(b));
    }

    static TrieNode* loadTrieNode(ifstream& in, TrieNode* parent, TrieIndex* index) {
        int marker;
        in.read(reinterpret_cast<char*>(&marker), sizeof(marker));
        if (marker == -1) return nullptr;

        int layer, entry, flag, size, other_size, actual_size, graph_size;
        bool has_index;
        in.read(reinterpret_cast<char*>(&layer), sizeof(layer));
        in.read(reinterpret_cast<char*>(&entry), sizeof(entry));
        in.read(reinterpret_cast<char*>(&flag), sizeof(flag));
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        in.read(reinterpret_cast<char*>(&other_size), sizeof(other_size));
        in.read(reinterpret_cast<char*>(&actual_size), sizeof(actual_size));
        in.read(reinterpret_cast<char*>(&graph_size), sizeof(graph_size));
        in.read(reinterpret_cast<char*>(&has_index), sizeof(has_index));

        TrieNode* node = new TrieNode(parent, layer);
        node->entry = entry;
        node->flag = flag;
        node->other_size = other_size;
        node->actual_size = actual_size;
        node->graph_size = graph_size;
        node->size = size;
        node->son = nullptr;
        if (has_index || node->parent == nullptr) {
            node->pre_index = node;
        } else {
            node->pre_index = node->parent->pre_index;
        }

        size_t element_count;
        in.read(reinterpret_cast<char*>(&element_count), sizeof(element_count));
        node->ids.resize(element_count);
        for (size_t i = 0; i < element_count; i++) {
            int id;
            in.read(reinterpret_cast<char*>(&id), sizeof(id));
            node->ids[i] = id;
        }

        size_t child_count;
        in.read(reinterpret_cast<char*>(&child_count), sizeof(child_count));
        for (size_t i = 0; i < child_count; i++) {
            int key;
            in.read(reinterpret_cast<char*>(&key), sizeof(key));
            node->child[key] = loadTrieNode(in, node, index);
            if(getlog(node->size,2) == getlog(node->child[key]->size,2))
                node->son = node->child[key];
        }
        return node;
    }


    __m128 masked_read(int dim, const float *x) const {
        __attribute__((__aligned__(16))) float buf[4] = {0, 0, 0, 0};
        switch (dim) {
            case 3:
                buf[2] = x[2];
            case 2:
                buf[1] = x[1];
            case 1:
                buf[0] = x[0];
        }
        return _mm_load_ps(buf);
    }

};

#endif //LABELFILTEREDANNS_TRIE_HPP
