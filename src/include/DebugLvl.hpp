#ifndef DEBUGLVL_H
#define DEBUGLVL_H

class DebugLevel {
   private:
    inline static int _debug_lvl{-1};

   public:
    DebugLevel() = default;
    ~DebugLevel() = default;
    DebugLevel(const DebugLevel&) = default;
    DebugLevel(DebugLevel&&) = default;
    auto        operator=(const DebugLevel&) -> DebugLevel& = default;
    auto        operator=(DebugLevel&&) -> DebugLevel& = default;

    friend void set_debug_lvl(int _lvl);
    friend auto debug_lvl(int _lvl) -> bool;
};

inline auto debug_lvl(int _lvl) -> bool {
    return (DebugLevel::_debug_lvl > _lvl);
}

inline void set_debug_lvl(int _lvl) {
    DebugLevel::_debug_lvl = _lvl;
}

#endif  // DEBUGLVL_H