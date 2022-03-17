// 分岐位置、分岐位置以降の葉の数、次の遷移先を保存
struct info_fp {
    uint64_t match; // 分岐位置
    uint64_t cnt; // 分岐位置以降の葉の数
    // uint64_t next_id; // 次のid
    std::vector<std::pair<uint64_t, uint64_t>> children; // 子の集合(左：next_id ，右：next_idの葉の数)

    info_fp() : match(0), cnt(0) {}

    info_fp(uint64_t m, uint64_t c) : match(m), cnt(c) {}
};

// 親情報を保存するときに使用
struct parent_info {
    uint64_t parent;
    uint64_t match;
    uint8_t c;
};