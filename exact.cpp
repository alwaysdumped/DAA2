#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <map>
#include <set>
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
    vector<vector<Edge>> g;
    vector<int> level, ptr;
    Dinic(int n) : g(n), level(n), ptr(n) {}
    void addEdge(int from, int to, int cap) {
        g[from].emplace_back(to, g[to].size(), cap);
        g[to].emplace_back(from, g[from].size() - 1, 0);
    }
    bool bfs(int s, int t) {
        fill(level.begin(), level.end(), -1);
        queue<int> q; q.push(s); level[s] = 0;
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (const Edge& e : g[v]) {
                if (e.cap > 0 && level[e.to] == -1) {
                    level[e.to] = level[v] + 1;
                    q.push(e.to);
                }
            }
        }
        return level[t] != -1;
    }
    int dfs(int v, int t, int flow) {
        if (v == t) return flow;
        for (int& i = ptr[v]; i < (int)g[v].size(); ++i) {
            Edge& e = g[v][i];
            if (e.cap > 0 && level[e.to] == level[v] + 1) {
                int pushed = dfs(e.to, t, min(flow, e.cap));
                if (pushed) {
                    e.cap -= pushed;
                    g[e.to][e.rev].cap += pushed;
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
            while (int pushed = dfs(s, t, INT_MAX)) flow += pushed;
        }
        return flow;
    }
    vector<bool> minCut(int s) {
        vector<bool> vis(g.size(), false);
        queue<int> q; q.push(s); vis[s] = true;
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (const Edge& e : g[v]) {
                if (e.cap > 0 && !vis[e.to]) {
                    vis[e.to] = true;
                    q.push(e.to);
                }
            }
        }
        return vis;
    }
};

struct Graph {
    int V;
    vector<vector<int>> adj;
    Graph(int V) : V(V), adj(V) {}
    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

// Motif enumeration
vector<vector<int>> enumerateEdges(const Graph& g) {
    vector<vector<int>> edges;
    for (int u = 0; u < g.V; ++u)
        for (int v : g.adj[u])
            if (u < v) edges.push_back({u, v});
    return edges;
}
vector<vector<int>> enumerateTriangles(const Graph& g) {
    vector<vector<int>> triangles;
    for (int u = 0; u < g.V; ++u)
        for (int v : g.adj[u]) if (v > u)
            for (int w : g.adj[v]) if (w > v && find(g.adj[u].begin(), g.adj[u].end(), w) != g.adj[u].end())
                triangles.push_back({u, v, w});
    return triangles;
}
void bronKerbosch(vector<vector<int>>& cliques, vector<int>& R, vector<int>& P, vector<int>& X, const vector<vector<int>>& adj, int h) {
    if (R.size() == h) { cliques.push_back(R); return; }
    if (P.empty() && X.empty()) return;
    vector<int> P_copy = P;
    for (int v : P_copy) {
        R.push_back(v);
        vector<int> newP, newX;
        for (int u : P)
            if (find(adj[v].begin(), adj[v].end(), u) != adj[v].end()) newP.push_back(u);
        for (int u : X)
            if (find(adj[v].begin(), adj[v].end(), u) != adj[v].end()) newX.push_back(u);
        bronKerbosch(cliques, R, newP, newX, adj, h);
        R.pop_back();
        P.erase(find(P.begin(), P.end(), v));
        X.push_back(v);
    }
}
vector<vector<int>> enumerateKCliques(const Graph& g, int h) {
    vector<vector<int>> cliques;
    vector<int> R, P, X;
    for (int i = 0; i < g.V; ++i) P.push_back(i);
    bronKerbosch(cliques, R, P, X, g.adj, h);
    return cliques;
}
vector<vector<int>> enumerate2Stars(const Graph& g) {
    vector<vector<int>> stars;
    for (int v = 0; v < g.V; ++v) {
        int d = g.adj[v].size();
        for (int i = 0; i < d; ++i)
            for (int j = i+1; j < d; ++j)
                stars.push_back({v, g.adj[v][i], g.adj[v][j]});
    }
    return stars;
}
vector<vector<int>> enumerateDiamonds(const Graph& g) {
    set<vector<int>> diamonds;
    for (int u = 0; u < g.V; ++u) {
        for (int v : g.adj[u]) if (v > u) {
            // Find all triangles on (u,v)
            vector<int> tri;
            for (int w : g.adj[u]) if (w != v && find(g.adj[v].begin(), g.adj[v].end(), w) != g.adj[v].end()) tri.push_back(w);
            int sz = tri.size();
            for (int i = 0; i < sz; ++i)
                for (int j = i+1; j < sz; ++j) {
                    vector<int> d = {u, v, tri[i], tri[j]};
                    sort(d.begin(), d.end());
                    diamonds.insert(d);
                }
        }
    }
    return vector<vector<int>>(diamonds.begin(), diamonds.end());
}

// Algorithm 1: Exact Densest Subgraph Discovery for a motif
pair<vector<int>, double> exactDSD(const Graph& g, const string& pattern) {
    int motif_size = 0;
    vector<vector<int>> motifs;
    if (pattern == "edge") { motif_size = 2; motifs = enumerateEdges(g); }
    else if (pattern == "triangle") { motif_size = 3; motifs = enumerateTriangles(g); }
    else if (pattern == "4clique") { motif_size = 4; motifs = enumerateKCliques(g, 4); }
    else if (pattern == "5clique") { motif_size = 5; motifs = enumerateKCliques(g, 5); }
    else if (pattern == "6clique") { motif_size = 6; motifs = enumerateKCliques(g, 6); }
    else if (pattern == "2star") { motif_size = -2; motifs = enumerate2Stars(g); }
    else if (pattern == "diamond") { motif_size = -3; motifs = enumerateDiamonds(g); }
    else { cout << "Unknown pattern: " << pattern << endl; exit(1); }

    int n = g.V;
    vector<int> motif_deg(n, 0);
    for (const auto& motif : motifs)
        for (int v : motif) motif_deg[v]++;
    double l = 0, u = *max_element(motif_deg.begin(), motif_deg.end());
    vector<int> best;
    double eps = 1.0 / (n * (n - 1) + 1);
    int max_iter = 30, iter = 0;
    while (u - l > eps && iter++ < max_iter) {
        double alpha = (l + u) / 2;
        int S = 0, T = 1;
        int N = 2 + n + motifs.size();
        Dinic dinic(N);
        for (int v = 0; v < n; ++v) dinic.addEdge(S, 2 + v, motif_deg[v]);
        for (int v = 0; v < n; ++v) dinic.addEdge(2 + v, T, int(alpha * abs(motif_size) + 0.5));
        for (int i = 0; i < motifs.size(); ++i) {
            int mnode = 2 + n + i;
            dinic.addEdge(mnode, T, 1);
            for (int v : motifs[i]) dinic.addEdge(2 + v, mnode, 1);
        }
        dinic.maxFlow(S, T);
        vector<bool> in_S = dinic.minCut(S);
        vector<int> sub;
        for (int v = 0; v < n; ++v)
            if (in_S[2 + v]) sub.push_back(v);
        if (!sub.empty()) {
            l = alpha;
            best = sub;
        } else {
            u = alpha;
        }
    }
    // Compute density for best subgraph
    double density = 0.0;
    if (!best.empty()) {
        map<int, int> o2n; for (int i = 0; i < best.size(); ++i) o2n[best[i]] = i;
        Graph finalg(best.size());
        for (int i = 0; i < best.size(); ++i) {
            int v = best[i];
            for (int u : g.adj[v])
                if (o2n.count(u) && o2n[u] > i)
                    finalg.addEdge(i, o2n[u]);
        }
        vector<vector<int>> motifs2;
        if (pattern == "edge") motifs2 = enumerateEdges(finalg);
        else if (pattern == "triangle") motifs2 = enumerateTriangles(finalg);
        else if (pattern == "4clique") motifs2 = enumerateKCliques(finalg, 4);
        else if (pattern == "5clique") motifs2 = enumerateKCliques(finalg, 5);
        else if (pattern == "6clique") motifs2 = enumerateKCliques(finalg, 6);
        else if (pattern == "2star") motifs2 = enumerate2Stars(finalg);
        else if (pattern == "diamond") motifs2 = enumerateDiamonds(finalg);
        density = double(motifs2.size()) / best.size();
    }
    return {best, density};
}

// File reading
Graph readGraphFromFile(const string& filename, int& n, int& m, bool& zero_based) {
    ifstream fin(filename.c_str());
    if (!fin) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    vector<pair<int, int>> edges;
    string line;
    int max_node = -1;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '%') continue;
        istringstream iss(line);
        int u, v;
        if (!(iss >> u >> v)) continue;
        edges.push_back(make_pair(u, v));
        max_node = max(max_node, max(u, v));
    }
    int min_node = INT_MAX;
    for (const auto& e : edges)
        min_node = min(min_node, min(e.first, e.second));
    zero_based = (min_node == 0);
    n = max_node + (zero_based ? 1 : 0);
    m = edges.size();
    Graph g(n);
    for (const auto& e : edges) {
        int u = e.first - (zero_based ? 0 : 1);
        int v = e.second - (zero_based ? 0 : 1);
        if (u >= 0 && u < n && v >= 0 && v < n)
            g.addEdge(u, v);
    }
    return g;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <graph_file> <pattern>" << endl;
        cout << "Patterns: edge, triangle, 4clique, 5clique, 6clique, 2star, diamond" << endl;
        return 1;
    }
    string filename = argv[1];
    string pattern = argv[2];
    int n, m;
    bool zero_based;
    clock_t start = clock();
    Graph g = readGraphFromFile(filename, n, m, zero_based);
    cout << "Graph has " << n << " nodes and " << m << " edges" << endl;
    cout << "Using " << (zero_based ? "0-based" : "1-based") << " indexing" << endl;
    cout << "Running Exact for pattern: " << pattern << "..." << endl;
    auto result = exactDSD(g, pattern);
    cout << "Exact Densest Subgraph (Density " << result.second << ", Size " << result.first.size() << "): ";
    for (int v : result.first) cout << (zero_based ? v : v+1) << " ";
    cout << endl;
    double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
    cout << "Total time: " << elapsed << " seconds" << endl;
    return 0;
}
