#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <map>
#include <fstream>
#include <sstream>
#include <climits>
#include <ctime>
using namespace std;

struct Edge {
    int to, rev, cap;
    Edge(int to, int rev, int cap) : to(to), rev(rev), cap(cap) {}
};

class Dinic {
public:
    vector<vector<Edge>> graph;
    vector<int> level, ptr;
    
    Dinic(int n) : graph(n), level(n), ptr(n) {}

    void addEdge(int from, int to, int cap) {
        if (from < 0 || to < 0 || from >= (int)graph.size() || to >= (int)graph.size()) return;
        graph[from].emplace_back(to, graph[to].size(), cap);
        graph[to].emplace_back(from, graph[from].size() - 1, 0);
    }

    bool bfs(int s, int t) {
        fill(level.begin(), level.end(), -1);
        queue<int> q;
        level[s] = 0;
        q.push(s);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (const Edge& e : graph[u]) {
                if (e.cap > 0 && level[e.to] == -1) {
                    level[e.to] = level[u] + 1;
                    q.push(e.to);
                }
            }
        }
        return level[t] != -1;
    }

    int dfs(int u, int t, int flow) {
        if (u == t) return flow;
        for (int& i = ptr[u]; i < (int)graph[u].size(); ++i) {
            Edge& e = graph[u][i];
            if (e.cap > 0 && level[e.to] == level[u] + 1) {
                int pushed = dfs(e.to, t, min(flow, e.cap));
                if (pushed > 0) {
                    e.cap -= pushed;
                    graph[e.to][e.rev].cap += pushed;
                    return pushed;
                }
            }
        }
        return 0;
    }

    int maxFlow(int s, int t) {
        int flow = 0;
        while (bfs(s, t)) {
            fill(ptr.begin(), ptr.end(), 0);
            while (int pushed = dfs(s, t, INT_MAX)) {
                flow += pushed;
            }
        }
        return flow;
    }
    
    vector<bool> minCut(int s) {
        vector<bool> in_S(graph.size(), false);
        queue<int> q;
        q.push(s);
        in_S[s] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (size_t i = 0; i < graph[u].size(); ++i) {
                const Edge& e = graph[u][i];
                if (e.cap > 0 && !in_S[e.to]) {
                    in_S[e.to] = true;
                    q.push(e.to);
                }
            }
        }
        return in_S;
    }
};

struct Graph {
    int V;
    vector<vector<int>> adj;
    
    Graph(int V) : V(V), adj(V) {}
    
    void addEdge(int u, int v) {
        if (u >= 0 && u < V && v >= 0 && v < V) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
    }
};

// Optimized triangle enumeration for h=3
vector<vector<int>> enumerateTriangles(const Graph& g) {
    vector<vector<int>> triangles;
    for (int u = 0; u < g.V; ++u) {
        for (int v : g.adj[u]) {
            if (v > u) { // Process each edge once
                for (int w : g.adj[v]) {
                    if (w > v) { // Ensure u < v < w
                        // Check if w is connected to u
                        if (find(g.adj[u].begin(), g.adj[u].end(), w) != g.adj[u].end()) {
                            triangles.push_back({u, v, w});
                        }
                    }
                }
            }
        }
    }
    return triangles;
}

// Optimized edge enumeration for h=2
vector<vector<int>> enumerateEdges(const Graph& g) {
    vector<vector<int>> edges;
    for (int u = 0; u < g.V; ++u) {
        for (int v : g.adj[u]) {
            if (v > u) { // Process each edge once
                edges.push_back({u, v});
            }
        }
    }
    return edges;
}

vector<vector<int>> enumerateCliques(const Graph& g, int h) {
    if (h == 2) return enumerateEdges(g);
    if (h == 3) return enumerateTriangles(g);
    
    // General case for other h values
    vector<vector<int>> cliques;
    vector<int> current;
    function<void(int, int)> backtrack = [&](int start, int remaining) {
        if (remaining == 0) {
            cliques.push_back(current);
            return;
        }
        for (int i = start; i < g.V; ++i) {
            bool canAdd = true;
            for (int u : current) {
                if (find(g.adj[u].begin(), g.adj[u].end(), i) == g.adj[u].end()) {
                    canAdd = false;
                    break;
                }
            }
            if (canAdd) {
                current.push_back(i);
                backtrack(i + 1, remaining - 1);
                current.pop_back();
            }
        }
    };
    backtrack(0, h);
    return cliques;
}

// Count clique degree for each vertex
vector<int> computeCliqueDegrees(const Graph& g, int h, const vector<vector<int>>& cliques) {
    vector<int> degree(g.V, 0);
    for (const auto& clique : cliques) {
        for (int v : clique) {
            degree[v]++;
        }
    }
    return degree;
}

// Algorithm 1: Exact
pair<vector<int>, double> exactDSD(const Graph& g, int h) {
    if (g.V <= 1) {
        vector<int> result;
        for (int i = 0; i < g.V; i++) result.push_back(i);
        return make_pair(result, 0.0);
    }
    
    int h_minus_1 = h - 1;
    cout << "Enumerating " << h_minus_1 << "-cliques for flow network..." << endl;
    auto A = enumerateCliques(g, h_minus_1);
    cout << "Found " << A.size() << " " << h_minus_1 << "-cliques" << endl;
    
    if (A.empty()) return make_pair(vector<int>(), 0.0);
    
    // Compute clique degrees
    vector<int> clique_degree(g.V, 0);
    for (int psi_idx = 0; psi_idx < (int)A.size(); ++psi_idx) {
        const vector<int>& psi = A[psi_idx];
        for (int v = 0; v < g.V; ++v) {
            bool is_clique = true;
            for (int u : psi) {
                if (u != v && find(g.adj[v].begin(), g.adj[v].end(), u) == g.adj[v].end()) {
                    is_clique = false;
                    break;
                }
            }
            if (is_clique) clique_degree[v]++;
        }
    }
    
    // Binary search bounds
    double low = 0.0, high = 0.0;
    for (int v = 0; v < g.V; ++v) {
        high = max(high, (double)clique_degree[v]);
    }
    
    cout << "Maximum clique degree: " << high << endl;
    
    if (high == 0) {
        return make_pair(vector<int>(), 0.0);
    }
    
    vector<int> best_subgraph;
    double epsilon = 1.0 / (g.V * (g.V - 1) + 1);
    
    cout << "Starting binary search with range [" << low << ", " << high << "]" << endl;
    
    while (high - low >= epsilon) {
        double alpha = (low + high) / 2;
        cout << "Trying alpha = " << alpha << endl;
        
        int s = 0, t = 1;
        int total_nodes = 2 + g.V + (int)A.size();
        Dinic dinic(total_nodes);
        
        // Source to vertices
        for (int v = 0; v < g.V; ++v) {
            dinic.addEdge(s, 2 + v, clique_degree[v]);
        }
        
        // Vertices to sink
        for (int v = 0; v < g.V; ++v) {
            dinic.addEdge(2 + v, t, (int)(alpha * h + 0.5));
        }
        
        // Clique nodes
        for (int psi_idx = 0; psi_idx < (int)A.size(); ++psi_idx) {
            int psi_node = 2 + g.V + psi_idx;
            dinic.addEdge(psi_node, t, 1);
            for (int v : A[psi_idx]) {
                if (v < g.V) {
                    dinic.addEdge(2 + v, psi_node, 1);
                }
            }
        }
        
        int flow = dinic.maxFlow(s, t);
        cout << "Max flow: " << flow << endl;
        
        // Find min-cut
        vector<bool> in_S = dinic.minCut(s);
        
        vector<int> subgraph;
        for (int v = 0; v < g.V; ++v) {
            if (2 + v < (int)in_S.size() && in_S[2 + v]) {
                subgraph.push_back(v);
            }
        }
        
        cout << "Subgraph size: " << subgraph.size() << " vertices" << endl;
        
        if (!subgraph.empty()) {
            low = alpha;
            best_subgraph = subgraph;
        } else {
            high = alpha;
        }
    }
    
    // Calculate actual density
    double density = 0.0;
    if (!best_subgraph.empty()) {
        // Create induced subgraph
        Graph sg(best_subgraph.size());
        map<int, int> old_to_new;
        for (size_t i = 0; i < best_subgraph.size(); ++i) {
            old_to_new[best_subgraph[i]] = i;
        }
        
        for (size_t i = 0; i < best_subgraph.size(); ++i) {
            int v = best_subgraph[i];
            for (int u : g.adj[v]) {
                auto it = old_to_new.find(u);
                if (it != old_to_new.end() && it->second > i) {
                    sg.addEdge(i, it->second);
                }
            }
        }
        
        auto sg_cliques = enumerateCliques(sg, h);
        density = (double)sg_cliques.size() / best_subgraph.size();
    }
    
    return make_pair(best_subgraph, density);
}

Graph readGraphFromFile(const string& filename, int& n, int& m, bool& zero_based) {
    ifstream fin(filename.c_str());
    if (!fin) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    
    vector<pair<int, int>> edges;
    string line;
    int max_node = -1;
    
    // Read all edges
    while (getline(fin, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#' || line[0] == '%') continue;
        
        istringstream iss(line);
        int u, v;
        if (!(iss >> u >> v)) continue; // Skip invalid lines
        
        edges.push_back(make_pair(u, v));
        max_node = max(max_node, max(u, v));
    }
    
    // Detect if 0-based or 1-based indexing
    int min_node = INT_MAX;
    for (const auto& e : edges) {
        min_node = min(min_node, min(e.first, e.second));
    }
    
    zero_based = (min_node == 0);
    n = max_node + (zero_based ? 1 : 0);
    m = edges.size();
    
    Graph g(n);
    for (const auto& e : edges) {
        int u = e.first - (zero_based ? 0 : 1);
        int v = e.second - (zero_based ? 0 : 1);
        if (u >= 0 && u < n && v >= 0 && v < n) {
            g.addEdge(u, v);
        }
    }
    
    return g;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <graph_file> <h>" << endl;
        return 1;
    }
    
    string filename = argv[1];
    int h = atoi(argv[2]);
    
    cout << "Reading graph from " << filename << "..." << endl;
    int n, m;
    bool zero_based;
    clock_t start = clock();
    
    try {
        Graph g = readGraphFromFile(filename, n, m, zero_based);
        cout << "Graph has " << n << " nodes and " << m << " edges" << endl;
        cout << "Using " << (zero_based ? "0-based" : "1-based") << " indexing" << endl;
        
        cout << "Computing densest subgraph with h=" << h << " using Algorithm 1 (Exact)..." << endl;
        pair<vector<int>, double> result = exactDSD(g, h);
        vector<int> densest_subgraph = result.first;
        double density = result.second;
        
        cout << "Exact Densest Subgraph (Density " << density << "): ";
        for (int v : densest_subgraph) {
            cout << (zero_based ? v : v+1) << " ";
        }
        cout << endl;
        
        double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
        cout << "Total time: " << elapsed << " seconds" << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
