#include "label_snapshot.hpp"
#include <iostream>
#include <stdexcept>

using namespace RouteOpt::Application::CVRP;

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    Label low, middle, high, inserted;
    low.rc = 2; middle.rc = 4; high.rc = 6; inserted.rc = 9;
    std::list<Label *> row0[2], row1[2];
    std::list<Label *> *buckets[]{row0, row1};
    row1[0] = {&low, &middle, &high};
    LabelSnapshot snapshot(buckets, 2, 2);
    row1[0].push_front(&inserted);
    row1[0].splice(row1[0].begin(), row1[0], std::prev(row1[0].end()));
    row1[0].clear();
    row1[1].push_back(&inserted);

    const auto saved = snapshot.bucket(1, 0);
    require(saved.size() == 3 && saved[0] == &low && saved[1] == &middle && saved[2] == &high,
            "Mutation of the live list changed the frozen partner order");
    require(snapshot.bucket(0, 0).empty() && snapshot.bucket(0, 1).empty() &&
            snapshot.bucket(1, 1).empty(), "Empty snapshot buckets changed after capture");
    LabelSnapshot current(buckets, 2, 2);
    require(current.bucket(1, 0).empty() && current.bucket(1, 1).front() == &inserted,
            "A new snapshot did not reflect the current buckets");
    row1[1].clear();
    LabelSnapshot empty(buckets, 2, 2);
    require(empty.labels.empty() && empty.bucket(1, 1).empty(), "All-empty snapshot failed");
    std::cout << "PASS: frozen order, live-list mutation, empty buckets, fresh capture\n";
}
