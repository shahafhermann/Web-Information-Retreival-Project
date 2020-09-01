package webdata.utils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 * A class for defining the probabilities for certain letters and errors.
 */
public class LetterProbability implements Serializable {

    //---------------------------------------------------//
    //***** Hard Coded Kernighan Confusion Matrices *****//
    //---------------------------------------------------//
    private static final String[] del_keys = {"tz", "tx", "ty", "tv", "tw", "tt", "tu", "tr", "ts", "tp", "tq", "tn", "to", "tl", "tm", "tj", "tk", "th", "ti", "tf", "tg", "td", "te", "tb", "tc", "ta", "me", "md", "mg", "mf", "ma", "mc", "mb", "mm", "ml", "mo", "mn", "mi", "mh", "mk", "mj", "mu", "mt", "mw", "mv", "mq", "mp", "ms", "mr", "my", "mx", "mz", "fp", "fq", "fr", "fs", "ft", "fu", "fv", "fw", "fx", "fy", "fz", "fa", "fb", "fc", "fd", "fe", "ff", "fg", "fh", "fi", "fj", "fk", "fl", "fm", "fn", "fo", "yi", "yh", "yk", "yj", "ym", "yl", "yo", "yn", "ya", "yc", "yb", "ye", "yd", "yg", "yf", "yy", "yx", "yz", "yq", "yp", "ys", "yr", "yu", "yt", "yw", "yv", "rt", "ru", "rv", "rw", "rp", "rq", "rr", "rs", "rx", "ry", "rz", "rd", "re", "rf", "rg", "ra", "rb", "rc", "rl", "rm", "rn", "ro", "rh", "ri", "rj", "rk", "kc", "kb", "ka", "kg", "kf", "ke", "kd", "kk", "kj", "ki", "kh", "ko", "kn", "km", "kl", "ks", "kr", "kq", "kp", "kw", "kv", "ku", "kt", "kz", "ky", "kx", "dn", "do", "dl", "dm", "dj", "dk", "dh", "di", "df", "dg", "dd", "de", "db", "dc", "da", "dz", "dx", "dy", "dv", "dw", "dt", "du", "dr", "ds", "dp", "dq", "wg", "wf", "we", "wd", "wc", "wb", "wa", "wo", "wn", "wm", "wl", "wk", "wj", "wi", "wh", "ww", "wv", "wu", "wt", "ws", "wr", "wq", "wp", "wz", "wy", "wx", "pr", "ps", "pp", "pq", "pv", "pw", "pt", "pu", "pz", "px", "py", "pb", "pc", "pa", "pf", "pg", "pd", "pe", "pj", "pk", "ph", "pi", "pn", "po", "pl", "pm", "iy", "ix", "iz", "iq", "ip", "is", "ir", "iu", "it", "iw", "iv", "ii", "ih", "ik", "ij", "im", "il", "io", "in", "ia", "ic", "ib", "ie", "id", "ig", "if", "bd", "be", "bf", "bg", "ba", "bb", "bc", "bl", "bm", "bn", "bo", "bh", "bi", "bj", "bk", "bt", "bu", "bv", "bw", "bp", "bq", "br", "bs", "bx", "by", "bz", "@g", "uy", "ux", "uz", "uu", "ut", "uw", "uv", "uq", "up", "us", "ur", "um", "ul", "uo", "un", "ui", "uh", "uk", "uj", "ue", "ud", "ug", "uf", "ua", "uc", "ub", "@e", "nh", "ni", "nj", "nk", "nl", "nm", "nn", "no", "na", "nb", "nc", "nd", "ne", "nf", "ng", "nx", "ny", "nz", "np", "nq", "nr", "ns", "nt", "nu", "nv", "nw", "gw", "gv", "gu", "gt", "gs", "gr", "gq", "gp", "gz", "gy", "gx", "gg", "gf", "ge", "gd", "gc", "gb", "ga", "go", "gn", "gm", "gl", "gk", "gj", "gi", "gh", "zl", "zm", "zn", "zo", "zh", "zi", "zj", "zk", "zd", "ze", "zf", "zg", "za", "zb", "zc", "zx", "zy", "zz", "zt", "zu", "zv", "zw", "zp", "zq", "zr", "zs", "@s", "@q", "sz", "sy", "sx", "ss", "sr", "sq", "sp", "sw", "sv", "su", "st", "sk", "sj", "si", "sh", "so", "sn", "sm", "sl", "sc", "sb", "sa", "sg", "sf", "se", "sd", "lf", "lg", "ld", "le", "lb", "lc", "la", "ln", "lo", "ll", "lm", "lj", "lk", "lh", "li", "lv", "lw", "lt", "lu", "lr", "ls", "lp", "lq", "lz", "lx", "ly", "em", "el", "eo", "en", "ei", "eh", "ek", "ej", "ee", "ed", "eg", "ef", "ea", "ec", "eb", "ey", "ex", "ez", "eu", "et", "ew", "ev", "eq", "ep", "es", "er", "xj", "xk", "xh", "xi", "xn", "xo", "xl", "xm", "xb", "xc", "xa", "xf", "xg", "xd", "xe", "xz", "xx", "xy", "xr", "xs", "xp", "xq", "xv", "xw", "xt", "xu", "@f", "qq", "qp", "qs", "qr", "qu", "qt", "qw", "qv", "qy", "qx", "qz", "qa", "qc", "qb", "qe", "qd", "qg", "qf", "qi", "qh", "qk", "qj", "qm", "ql", "qo", "qn", "@o", "@n", "@m", "@c", "@b", "@a", "@k", "@j", "@i", "@h", "jx", "jy", "jz", "@l", "jt", "ju", "jv", "jw", "jp", "jq", "jr", "js", "jl", "jm", "jn", "jo", "jh", "ji", "jj", "jk", "jd", "je", "jf", "jg", "@w", "ja", "jb", "jc", "@z", "@y", "@x", "ck", "cj", "ci", "ch", "co", "cn", "cm", "cl", "cc", "cb", "ca", "cg", "cf", "ce", "cd", "cz", "cy", "cx", "@r", "cs", "cr", "cq", "cp", "cw", "cv", "cu", "ct", "@d", "@p", "@v", "@u", "@t", "va", "vb", "vc", "vd", "ve", "vf", "vg", "vh", "vi", "vj", "vk", "vl", "vm", "vn", "vo", "vp", "vq", "vr", "vs", "vt", "vu", "vv", "vw", "vx", "vy", "vz", "oo", "on", "om", "ol", "ok", "oj", "oi", "oh", "og", "of", "oe", "od", "oc", "ob", "oa", "oz", "oy", "ox", "ow", "ov", "ou", "ot", "os", "or", "oq", "op", "hz", "hx", "hy", "hr", "hs", "hp", "hq", "hv", "hw", "ht", "hu", "hj", "hk", "hh", "hi", "hn", "ho", "hl", "hm", "hb", "hc", "ha", "hf", "hg", "hd", "he", "aa", "ac", "ab", "ae", "ad", "ag", "af", "ai", "ah", "ak", "aj", "am", "al", "ao", "an", "aq", "ap", "as", "ar", "au", "at", "aw", "av", "ay", "ax", "az"};
    private static final int[] del_vals = {0, 0, 2, 0, 4, 137, 14, 203, 5, 1, 0, 3, 11, 31, 3, 0, 0, 49, 427, 1, 7, 0, 76, 1, 2, 24, 33, 0, 0, 0, 15, 0, 10, 180, 0, 7, 7, 42, 1, 0, 0, 4, 0, 0, 0, 0, 31, 9, 0, 0, 0, 0, 0, 0, 11, 0, 8, 1, 0, 0, 0, 1, 0, 4, 0, 0, 0, 13, 46, 0, 0, 79, 0, 0, 12, 0, 0, 4, 1, 0, 0, 0, 2, 1, 1, 1, 2, 34, 1, 2, 0, 1, 0, 0, 0, 0, 0, 1, 17, 0, 0, 1, 1, 0, 68, 0, 10, 1, 2, 0, 277, 103, 0, 27, 0, 19, 188, 0, 11, 63, 4, 12, 33, 7, 157, 21, 5, 132, 0, 3, 0, 0, 4, 8, 1, 15, 1, 1, 0, 5, 1, 0, 17, 0, 3, 5, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3, 3, 8, 4, 1, 1, 0, 62, 0, 10, 25, 45, 0, 7, 12, 0, 0, 6, 2, 0, 0, 3, 11, 1, 0, 0, 0, 1, 11, 1, 0, 0, 40, 2, 2, 0, 1, 0, 0, 15, 11, 0, 0, 0, 0, 24, 2, 0, 0, 0, 0, 0, 58, 1, 93, 0, 0, 0, 18, 2, 0, 0, 0, 0, 0, 25, 0, 0, 0, 22, 0, 0, 12, 15, 0, 30, 28, 1, 1, 0, 7, 0, 7, 71, 16, 1, 64, 0, 1, 1, 0, 0, 0, 14, 38, 41, 82, 26, 60, 1, 23, 26, 9, 1, 0, 22, 0, 0, 2, 2, 1, 26, 0, 0, 2, 0, 183, 0, 0, 0, 6, 1, 0, 0, 0, 6, 17, 0, 0, 0, 7, 1, 0, 0, 0, 66, 0, 0, 0, 0, 31, 129, 2, 39, 1, 111, 28, 0, 0, 0, 15, 10, 1, 0, 26, 9, 6, 20, 0, 191, 0, 0, 0, 17, 144, 21, 21, 0, 42, 71, 68, 1, 160, 0, 2, 0, 0, 0, 0, 127, 87, 43, 1, 1, 0, 0, 22, 1, 7, 52, 0, 0, 0, 1, 0, 37, 1, 83, 2, 0, 0, 25, 4, 29, 0, 3, 0, 0, 39, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 26, 0, 0, 1, 0, 265, 4, 0, 30, 0, 0, 21, 124, 0, 0, 231, 18, 30, 0, 1, 2, 27, 0, 16, 0, 1, 74, 0, 0, 0, 6, 48, 0, 1, 24, 0, 29, 211, 2, 0, 0, 0, 217, 2, 0, 7, 3, 2, 12, 0, 0, 0, 0, 11, 9, 32, 19, 76, 6, 1, 0, 0, 89, 74, 1, 3, 80, 50, 1, 1, 7, 0, 8, 34, 1, 2, 1, 9, 223, 237, 0, 0, 1, 0, 0, 0, 0, 0, 0, 17, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 6, 0, 0, 0, 5, 0, 20, 0, 0, 0, 0, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 16, 41, 14, 20, 6, 3, 20, 6, 0, 0, 0, 22, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 24, 0, 0, 0, 2, 0, 0, 9, 0, 320, 24, 33, 0, 0, 17, 70, 0, 37, 0, 0, 63, 0, 0, 1, 0, 28, 6, 46, 0, 0, 0, 0, 17, 54, 31, 17, 1, 2, 6, 9, 0, 0, 0, 58, 0, 0, 0, 31, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 26, 70, 9, 13, 0, 1, 4, 0, 5, 0, 8, 6, 3, 4, 11, 0, 1, 0, 5, 2, 47, 13, 20, 98, 0, 20, 0, 0, 1, 15, 1, 0, 0, 0, 1, 26, 0, 0, 0, 25, 24, 9, 22, 7, 1, 12, 1, 15, 0, 0, 3, 20, 0, 58, 7, 3, 21, 18, 5, 61, 8, 4, 0, 5, 43, 0, 53, 0, 9, 28, 98, 62, 53, 0, 1, 2, 0, 0};

    private static final String[] ins_keys = {"gw", "gv", "gu", "gt", "gs", "gr", "gq", "gp", "@y", "gz", "gy", "gx", "gg", "gf", "ge", "gd", "gc", "gb", "ga", "go", "gn", "gm", "gl", "gk", "gj", "gi", "gh", "tz", "tx", "ty", "tv", "tw", "tt", "tu", "tr", "ts", "tp", "tq", "tn", "to", "tl", "tm", "tj", "tk", "th", "ti", "tf", "tg", "td", "te", "tb", "tc", "ta", "vu", "zl", "zm", "zn", "zo", "zh", "zi", "zj", "zk", "zd", "ze", "zf", "zg", "za", "zb", "zc", "zx", "zy", "zz", "zt", "zu", "zv", "zw", "zp", "zq", "zr", "zs", "@s", "wl", "@q", "va", "@c", "vc", "wk", "@p", "vh", "wj", "vi", "vj", "vk", "vl", "vm", "wi", "@v", "vn", "vo", "me", "md", "mg", "mf", "ma", "mc", "mb", "mm", "ml", "mo", "mn", "mi", "mh", "mk", "mj", "mu", "mt", "mw", "mv", "mq", "mp", "ms", "mr", "vt", "my", "mx", "mz", "vv", "vw", "@t", "vx", "vz", "@b", "fp", "fq", "fr", "fs", "ft", "fu", "fv", "fw", "fx", "fy", "fz", "fa", "fb", "fc", "fd", "fe", "ff", "fg", "fh", "fi", "fj", "fk", "fl", "fm", "fn", "fo", "sz", "sy", "sx", "ss", "sr", "sq", "sp", "sw", "sv", "su", "st", "sk", "sj", "si", "sh", "so", "sn", "sm", "sl", "sc", "sb", "sa", "sg", "sf", "se", "sd", "lf", "lg", "ld", "le", "lb", "lc", "la", "ln", "lo", "ll", "lm", "lj", "lk", "lh", "li", "lv", "lw", "lt", "lu", "lr", "ls", "lp", "lq", "lz", "lx", "ly", "wq", "yh", "yk", "yj", "ym", "yl", "yo", "yn", "ya", "yc", "yb", "ye", "yd", "yg", "yf", "yy", "yx", "yz", "yq", "yp", "ys", "yr", "yu", "yt", "yw", "yv", "em", "el", "eo", "en", "ei", "eh", "ek", "ej", "ee", "ed", "eg", "ef", "ea", "ec", "eb", "ey", "ex", "@g", "ez", "eu", "et", "ew", "ev", "eq", "ep", "es", "er", "rt", "ru", "rv", "rw", "rp", "rq", "rr", "rs", "rx", "ry", "rz", "rd", "re", "rf", "rg", "ra", "rb", "rc", "rl", "rm", "rn", "ro", "rh", "ri", "rj", "rk", "xj", "xk", "xh", "xi", "xn", "xo", "xl", "xm", "xb", "xc", "xa", "xf", "xg", "xd", "xe", "xz", "xx", "xy", "xr", "xs", "xp", "xq", "xv", "xw", "xt", "xu", "wy", "wx", "@d", "kc", "kb", "ka", "kg", "kf", "ke", "kd", "kk", "kj", "ki", "kh", "ko", "kn", "km", "kl", "ks", "kr", "kq", "kp", "kw", "kv", "ku", "kt", "kz", "ky", "kx", "dn", "do", "dl", "dm", "dj", "dk", "dh", "di", "df", "dg", "dd", "de", "db", "dc", "da", "dz", "dx", "dy", "dv", "dw", "dt", "du", "dr", "ds", "dp", "dq", "qq", "qp", "qs", "qr", "qu", "qt", "qw", "qv", "qy", "qx", "qz", "qa", "qc", "qb", "qe", "qd", "qg", "qf", "qi", "qh", "qk", "qj", "qm", "ql", "qo", "qn", "@k", "@j", "@e", "@i", "@h", "wc", "wb", "wa", "wo", "wn", "wm", "wg", "wf", "we", "wd", "jx", "jy", "jz", "@l", "jt", "ju", "jv", "jw", "jp", "jq", "jr", "js", "jl", "jm", "jn", "jo", "jh", "ji", "jj", "jk", "jd", "je", "jf", "jg", "@w", "ja", "jb", "jc", "ww", "wv", "wu", "wt", "ws", "wr", "ck", "cj", "ci", "ch", "co", "cn", "cm", "cl", "cc", "cb", "ca", "wp", "cg", "cf", "ce", "cd", "cz", "cy", "cx", "@r", "cs", "cr", "cq", "cp", "cw", "cv", "cu", "ct", "pr", "ps", "pp", "pq", "pv", "pw", "pt", "pu", "pz", "px", "py", "wz", "pb", "pc", "pa", "pf", "pg", "pd", "pe", "pj", "pk", "ph", "pi", "pn", "po", "pl", "pm", "iy", "ix", "vb", "iz", "vd", "ve", "vf", "vg", "iq", "ip", "is", "ir", "iu", "it", "iw", "iv", "ii", "ih", "ik", "ij", "im", "il", "io", "in", "ia", "vy", "ic", "ib", "ie", "id", "ig", "if", "@x", "wh", "yi", "@u", "vr", "@f", "@o", "@n", "@m", "vs", "bd", "be", "bf", "bg", "ba", "bb", "bc", "bl", "bm", "bn", "bo", "bh", "bi", "bj", "bk", "bt", "bu", "bv", "bw", "bp", "bq", "br", "bs", "bx", "by", "bz", "oo", "on", "om", "ol", "ok", "oj", "oi", "oh", "og", "of", "oe", "od", "oc", "ob", "oa", "oz", "oy", "ox", "ow", "ov", "ou", "ot", "os", "or", "oq", "op", "@a", "hz", "hx", "hy", "hr", "hs", "hp", "hq", "hv", "hw", "ht", "hu", "hj", "hk", "hh", "hi", "hn", "ho", "hl", "hm", "hb", "hc", "ha", "hf", "hg", "hd", "he", "uy", "ux", "uz", "uu", "ut", "uw", "uv", "uq", "up", "us", "ur", "um", "ul", "uo", "un", "ui", "uh", "uk", "uj", "ue", "ud", "ug", "uf", "ua", "uc", "ub", "aa", "ac", "ab", "ae", "ad", "ag", "af", "ai", "ah", "ak", "aj", "am", "al", "ao", "an", "aq", "ap", "as", "ar", "au", "at", "aw", "av", "ay", "ax", "az", "nh", "ni", "nj", "nk", "nl", "nm", "nn", "no", "na", "nb", "nc", "nd", "ne", "nf", "ng", "nx", "ny", "nz", "np", "nq", "nr", "ns", "nt", "nu", "nv", "nw", "vp", "@z", "vq"};
    private static final int[] ins_vals = {1, 0, 3, 2, 69, 5, 1, 0, 1, 0, 0, 0, 5, 1, 5, 0, 0, 0, 8, 1, 1, 0, 2, 0, 0, 8, 12, 0, 0, 6, 0, 5, 183, 11, 54, 264, 1, 0, 1, 23, 6, 3, 1, 0, 24, 59, 1, 10, 3, 65, 0, 0, 39, 1, 0, 0, 0, 0, 0, 6, 0, 0, 0, 5, 1, 0, 2, 0, 0, 0, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 23, 1, 0, 0, 9, 0, 1, 10, 0, 0, 10, 0, 0, 1, 0, 1, 1, 1, 0, 17, 0, 0, 0, 11, 1, 1, 102, 0, 7, 44, 6, 1, 1, 0, 2, 1, 1, 0, 0, 2, 47, 0, 0, 0, 0, 0, 5, 1, 2, 0, 0, 8, 0, 0, 5, 23, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2, 27, 1, 0, 12, 0, 0, 10, 0, 0, 0, 0, 7, 0, 205, 1, 0, 1, 1, 0, 7, 49, 2, 0, 101, 50, 3, 7, 10, 2, 7, 1, 13, 1, 0, 41, 20, 0, 0, 1, 38, 1, 0, 3, 0, 7, 128, 1, 0, 2, 0, 79, 1, 0, 7, 3, 0, 97, 0, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 6, 5, 2, 1, 3, 0, 0, 0, 2, 0, 0, 0, 0, 33, 1, 13, 1, 1, 0, 6, 4, 5, 27, 4, 1, 3, 0, 147, 76, 0, 2, 39, 8, 2, 8, 2, 14, 0, 4, 6, 10, 1, 0, 1, 417, 83, 29, 7, 0, 1, 0, 0, 132, 273, 0, 10, 0, 0, 89, 1, 1, 15, 2, 1, 5, 9, 7, 10, 2, 64, 0, 0, 0, 0, 6, 1, 0, 3, 0, 1, 0, 18, 0, 0, 0, 0, 1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 4, 2, 0, 0, 9, 1, 1, 0, 1, 1, 2, 0, 0, 1, 95, 0, 0, 1, 0, 0, 1, 0, 0, 4, 0, 9, 13, 6, 1, 0, 0, 0, 9, 2, 0, 17, 14, 0, 3, 18, 0, 0, 5, 0, 0, 0, 0, 6, 119, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 17, 1, 26, 5, 3, 0, 0, 0, 0, 2, 0, 0, 0, 10, 1, 0, 0, 0, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 4, 0, 2, 0, 8, 1, 3, 0, 50, 18, 7, 1, 1, 1, 54, 0, 19, 0, 0, 0, 13, 1, 0, 0, 1, 6, 25, 7, 0, 1, 0, 4, 8, 7, 29, 52, 70, 0, 1, 1, 9, 1, 0, 0, 0, 0, 0, 1, 23, 0, 0, 1, 10, 0, 0, 20, 3, 0, 26, 2, 0, 0, 1, 2, 1, 0, 36, 0, 0, 0, 1, 30, 9, 11, 29, 0, 0, 69, 1, 1, 2, 11, 17, 27, 33, 10, 0, 13, 3, 25, 13, 1, 0, 1, 1, 2, 11, 0, 11, 2, 2, 6, 0, 0, 7, 0, 1, 3, 11, 0, 15, 0, 1, 1, 0, 50, 0, 0, 0, 0, 3, 0, 0, 0, 5, 16, 0, 0, 0, 64, 13, 3, 6, 0, 1, 28, 0, 1, 2, 7, 3, 1, 1, 14, 1, 1, 0, 0, 1, 19, 4, 59, 16, 0, 30, 46, 0, 0, 3, 16, 24, 0, 0, 0, 5, 22, 1, 2, 0, 18, 17, 1, 4, 1, 0, 1, 0, 4, 0, 10, 1, 24, 3, 2, 0, 26, 27, 0, 0, 0, 3, 19, 49, 3, 3, 1, 9, 24, 1, 1, 1, 9, 0, 0, 0, 15, 3, 0, 15, 14, 1, 10, 7, 1, 0, 33, 1, 4, 1, 2, 31, 12, 39, 3, 4, 134, 28, 28, 7, 1, 0, 4, 1, 1, 0, 34, 0, 1, 1, 26, 99, 12, 15, 5, 7, 13, 52, 4, 17, 0, 1, 0, 0, 0, 2, 156, 53, 1, 1, 0, 1, 2, 0};

    private static final String[] sub_keys = {"gw", "gv", "gu", "gt", "gs", "gr", "gq", "gp", "gz", "gy", "gx", "gg", "gf", "ge", "gd", "gc", "gb", "ga", "go", "gn", "gm", "gl", "gk", "gj", "gi", "gh", "tz", "tx", "ty", "tv", "tw", "tt", "tu", "tr", "ts", "tp", "tq", "tn", "to", "tl", "tm", "tj", "tk", "th", "ti", "tf", "tg", "td", "te", "tb", "tc", "ta", "vu", "zl", "zm", "zn", "zo", "zh", "zi", "zj", "zk", "zd", "ze", "zf", "zg", "za", "zb", "zc", "zx", "zy", "zz", "zt", "zu", "zv", "zw", "zp", "zq", "zr", "zs", "wl", "va", "vc", "wk", "vh", "wj", "vi", "vj", "vk", "vl", "vm", "wi", "vn", "vo", "me", "md", "mg", "mf", "ma", "mc", "mb", "mm", "ml", "mo", "mn", "mi", "mh", "mk", "mj", "mu", "mt", "mw", "mv", "mq", "mp", "ms", "mr", "vt", "my", "mx", "mz", "vv", "vw", "vx", "vz", "fp", "fq", "fr", "fs", "ft", "fu", "fv", "fw", "fx", "fy", "fz", "fa", "fb", "fc", "fd", "fe", "ff", "fg", "fh", "fi", "fj", "fk", "fl", "fm", "fn", "fo", "sz", "sy", "sx", "ss", "sr", "sq", "sp", "sw", "sv", "su", "st", "sk", "sj", "si", "sh", "so", "sn", "sm", "sl", "sc", "sb", "sa", "sg", "sf", "se", "sd", "lf", "lg", "ld", "le", "lb", "lc", "la", "ln", "lo", "ll", "lm", "lj", "lk", "lh", "li", "lv", "lw", "lt", "lu", "lr", "ls", "lp", "lq", "lz", "lx", "ly", "wq", "yh", "yk", "yj", "ym", "yl", "yo", "yn", "ya", "yc", "yb", "ye", "yd", "yg", "yf", "yy", "yx", "yz", "yq", "yp", "ys", "yr", "yu", "yt", "yw", "yv", "em", "el", "eo", "en", "ei", "eh", "ek", "ej", "ee", "ed", "eg", "ef", "ea", "ec", "eb", "ey", "ex", "ez", "eu", "et", "ew", "ev", "eq", "ep", "es", "er", "rt", "ru", "rv", "rw", "rp", "rq", "rr", "rs", "rx", "ry", "rz", "rd", "re", "rf", "rg", "ra", "rb", "rc", "rl", "rm", "rn", "ro", "rh", "ri", "rj", "rk", "xj", "xk", "xh", "xi", "xn", "xo", "xl", "xm", "xb", "xc", "xa", "xf", "xg", "xd", "xe", "xz", "xx", "xy", "xr", "xs", "xp", "xq", "xv", "xw", "xt", "xu", "wy", "wx", "kc", "kb", "ka", "kg", "kf", "ke", "kd", "kk", "kj", "ki", "kh", "ko", "kn", "km", "kl", "ks", "kr", "kq", "kp", "kw", "kv", "ku", "kt", "kz", "ky", "kx", "dn", "do", "dl", "dm", "dj", "dk", "dh", "di", "df", "dg", "dd", "de", "db", "dc", "da", "dz", "dx", "dy", "dv", "dw", "dt", "du", "dr", "ds", "dp", "dq", "qq", "qp", "qs", "qr", "qu", "qt", "qw", "qv", "qy", "qx", "qz", "qa", "qc", "qb", "qe", "qd", "qg", "qf", "qi", "qh", "qk", "qj", "qm", "ql", "qo", "qn", "wc", "wb", "wa", "wo", "wn", "wm", "wg", "wf", "we", "wd", "jx", "jy", "jz", "jt", "ju", "jv", "jw", "jp", "jq", "jr", "js", "jl", "jm", "jn", "jo", "jh", "ji", "jj", "jk", "jd", "je", "jf", "jg", "ja", "jb", "jc", "ww", "wv", "wu", "wt", "ws", "wr", "ck", "cj", "ci", "ch", "co", "cn", "cm", "cl", "cc", "cb", "ca", "wp", "cg", "cf", "ce", "cd", "cz", "cy", "cx", "cs", "cr", "cq", "cp", "cw", "cv", "cu", "ct", "pr", "ps", "pp", "pq", "pv", "pw", "pt", "pu", "pz", "px", "py", "wz", "pb", "pc", "pa", "pf", "pg", "pd", "pe", "pj", "pk", "ph", "pi", "pn", "po", "pl", "pm", "iy", "ix", "vb", "iz", "vd", "ve", "vf", "vg", "iq", "ip", "is", "ir", "iu", "it", "iw", "iv", "ii", "ih", "ik", "ij", "im", "il", "io", "in", "ia", "vy", "ic", "ib", "ie", "id", "ig", "if", "wh", "yi", "vr", "vs", "bd", "be", "bf", "bg", "ba", "bb", "bc", "bl", "bm", "bn", "bo", "bh", "bi", "bj", "bk", "bt", "bu", "bv", "bw", "bp", "bq", "br", "bs", "bx", "by", "bz", "oo", "on", "om", "ol", "ok", "oj", "oi", "oh", "og", "of", "oe", "od", "oc", "ob", "oa", "oz", "oy", "ox", "ow", "ov", "ou", "ot", "os", "or", "oq", "op", "hz", "hx", "hy", "hr", "hs", "hp", "hq", "hv", "hw", "ht", "hu", "hj", "hk", "hh", "hi", "hn", "ho", "hl", "hm", "hb", "hc", "ha", "hf", "hg", "hd", "he", "uy", "ux", "uz", "uu", "ut", "uw", "uv", "uq", "up", "us", "ur", "um", "ul", "uo", "un", "ui", "uh", "uk", "uj", "ue", "ud", "ug", "uf", "ua", "uc", "ub", "aa", "ac", "ab", "ae", "ad", "ag", "af", "ai", "ah", "ak", "aj", "am", "al", "ao", "an", "aq", "ap", "as", "ar", "au", "at", "aw", "av", "ay", "ax", "az", "nh", "ni", "nj", "nk", "nl", "nm", "nn", "no", "na", "nb", "nc", "nd", "ne", "nf", "ng", "nx", "ny", "nz", "np", "nq", "nr", "ns", "nt", "nu", "nv", "nw", "vp", "vq"};
    private static final int[] sub_vals = {1, 0, 0, 21, 13, 5, 3, 1, 0, 3, 0, 0, 2, 9, 11, 11, 1, 4, 2, 0, 0, 3, 1, 1, 0, 0, 6, 0, 7, 2, 19, 0, 0, 11, 37, 6, 0, 5, 5, 14, 9, 1, 0, 5, 0, 5, 19, 42, 7, 4, 9, 3, 0, 7, 5, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 3, 0, 3, 0, 0, 0, 0, 0, 2, 21, 0, 0, 7, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 8, 0, 2, 1, 7, 3, 0, 4, 0, 180, 0, 6, 4, 0, 13, 15, 2, 3, 0, 6, 9, 0, 3, 3, 2, 0, 0, 0, 0, 0, 0, 0, 6, 4, 12, 0, 0, 2, 0, 0, 0, 0, 15, 0, 3, 1, 0, 5, 2, 0, 0, 0, 3, 4, 1, 0, 1, 20, 3, 0, 14, 0, 7, 5, 0, 0, 15, 0, 1, 0, 1, 1, 6, 0, 27, 27, 8, 11, 0, 4, 35, 33, 4, 5, 4, 0, 10, 1, 2, 14, 2, 0, 0, 0, 1, 6, 13, 0, 0, 2, 0, 11, 10, 5, 0, 0, 0, 0, 0, 7, 0, 0, 2, 0, 6, 0, 0, 2, 0, 15, 0, 1, 0, 0, 1, 0, 0, 1, 36, 7, 5, 8, 0, 0, 0, 3, 93, 5, 89, 0, 0, 0, 0, 11, 2, 2, 388, 3, 0, 18, 0, 0, 15, 6, 1, 0, 0, 0, 12, 14, 22, 4, 0, 0, 14, 0, 0, 12, 1, 0, 0, 30, 12, 2, 2, 0, 14, 0, 8, 4, 20, 1, 8, 2, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 8, 2, 1, 2, 1, 1, 4, 0, 0, 0, 5, 2, 0, 5, 0, 6, 0, 0, 0, 4, 0, 0, 0, 3, 0, 0, 3, 0, 3, 7, 0, 2, 5, 0, 0, 5, 0, 12, 10, 13, 1, 0, 0, 2, 0, 4, 22, 0, 43, 30, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 2, 1, 0, 0, 0, 0, 0, 0, 9, 0, 0, 1, 0, 1, 1, 0, 0, 1, 3, 3, 6, 1, 0, 0, 0, 1, 9, 7, 0, 0, 5, 6, 7, 5, 9, 0, 16, 0, 1, 1, 39, 5, 2, 10, 7, 3, 1, 40, 1, 3, 0, 0, 4, 1, 6, 0, 0, 0, 0, 0, 11, 1, 0, 6, 5, 2, 0, 9, 0, 0, 2, 6, 15, 2, 7, 15, 1, 0, 0, 0, 0, 3, 0, 0, 0, 2, 0, 47, 1, 2, 0, 0, 0, 0, 0, 0, 6, 49, 0, 103, 0, 0, 0, 146, 0, 1, 0, 2, 15, 0, 8, 9, 2, 2, 3, 0, 0, 9, 5, 11, 5, 0, 1, 0, 0, 0, 1, 0, 0, 8, 10, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 25, 0, 0, 0, 116, 3, 1, 1, 91, 0, 18, 0, 0, 0, 39, 14, 4, 2, 0, 14, 0, 0, 0, 3, 1, 3, 0, 0, 2, 11, 0, 0, 2, 0, 0, 14, 2, 0, 12, 8, 0, 1, 0, 0, 3, 0, 8, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 0, 0, 43, 2, 64, 0, 0, 0, 44, 0, 0, 0, 20, 0, 0, 0, 7, 0, 342, 1, 0, 0, 118, 2, 1, 0, 0, 0, 76, 3, 0, 0, 35, 1, 9, 9, 1, 0, 5, 0, 0, 19, 1, 0, 4, 35, 78, 0, 0, 2, 7, 6, 5, 3, 0, 1, 2, 0, 2, 7, 0, 28, 5, 7, 0, 0, 1, 0, 0};

    private static final String[] trans_keys = {"gw", "gv", "gu", "gt", "gs", "gr", "gq", "gp", "gz", "gy", "gx", "gg", "gf", "ge", "gd", "gc", "gb", "ga", "go", "gn", "gm", "gl", "gk", "gj", "gi", "gh", "tz", "tx", "ty", "tv", "tw", "tt", "tu", "tr", "ts", "tp", "tq", "tn", "to", "tl", "tm", "tj", "tk", "th", "ti", "tf", "tg", "td", "te", "tb", "tc", "ta", "vu", "zl", "zm", "zn", "zo", "zh", "zi", "zj", "zk", "zd", "ze", "zf", "zg", "za", "zb", "zc", "zx", "zy", "zz", "zt", "zu", "zv", "zw", "zp", "zq", "zr", "zs", "wl", "va", "vc", "wk", "vh", "wj", "vi", "vj", "vk", "vl", "vm", "wi", "vn", "vo", "me", "md", "mg", "mf", "ma", "mc", "mb", "mm", "ml", "mo", "mn", "mi", "mh", "mk", "mj", "mu", "mt", "mw", "mv", "mq", "mp", "ms", "mr", "vt", "my", "mx", "mz", "vv", "vw", "vx", "vz", "fp", "fq", "fr", "fs", "ft", "fu", "fv", "fw", "fx", "fy", "fz", "fa", "fb", "fc", "fd", "fe", "ff", "fg", "fh", "fi", "fj", "fk", "fl", "fm", "fn", "fo", "sz", "sy", "sx", "ss", "sr", "sq", "sp", "sw", "sv", "su", "st", "sk", "sj", "si", "sh", "so", "sn", "sm", "sl", "sc", "sb", "sa", "sg", "sf", "se", "sd", "lf", "lg", "ld", "le", "lb", "lc", "la", "ln", "lo", "ll", "lm", "lj", "lk", "lh", "li", "lv", "lw", "lt", "lu", "lr", "ls", "lp", "lq", "lz", "lx", "ly", "wq", "yh", "yk", "yj", "ym", "yl", "yo", "yn", "ya", "yc", "yb", "ye", "yd", "yg", "yf", "yy", "yx", "yz", "yq", "yp", "ys", "yr", "yu", "yt", "yw", "yv", "em", "el", "eo", "en", "ei", "eh", "ek", "ej", "ee", "ed", "eg", "ef", "ea", "ec", "eb", "ey", "ex", "ez", "eu", "et", "ew", "ev", "eq", "ep", "es", "er", "rt", "ru", "rv", "rw", "rp", "rq", "rr", "rs", "rx", "ry", "rz", "rd", "re", "rf", "rg", "ra", "rb", "rc", "rl", "rm", "rn", "ro", "rh", "ri", "rj", "rk", "xj", "xk", "xh", "xi", "xn", "xo", "xl", "xm", "xb", "xc", "xa", "xf", "xg", "xd", "xe", "xz", "xx", "xy", "xr", "xs", "xp", "xq", "xv", "xw", "xt", "xu", "wy", "wx", "kc", "kb", "ka", "kg", "kf", "ke", "kd", "kk", "kj", "ki", "kh", "ko", "kn", "km", "kl", "ks", "kr", "kq", "kp", "kw", "kv", "ku", "kt", "kz", "ky", "kx", "dn", "do", "dl", "dm", "dj", "dk", "dh", "di", "df", "dg", "dd", "de", "db", "dc", "da", "dz", "dx", "dy", "dv", "dw", "dt", "du", "dr", "ds", "dp", "dq", "qq", "qp", "qs", "qr", "qu", "qt", "qw", "qv", "qy", "qx", "qz", "qa", "qc", "qb", "qe", "qd", "qg", "qf", "qi", "qh", "qk", "qj", "qm", "ql", "qo", "qn", "wc", "wb", "wa", "wo", "wn", "wm", "wg", "wf", "we", "wd", "jx", "jy", "jz", "jt", "ju", "jv", "jw", "jp", "jq", "jr", "js", "jl", "jm", "jn", "jo", "jh", "ji", "jj", "jk", "jd", "je", "jf", "jg", "ja", "jb", "jc", "ww", "wv", "wu", "wt", "ws", "wr", "ck", "cj", "ci", "ch", "co", "cn", "cm", "cl", "cc", "cb", "ca", "wp", "cg", "cf", "ce", "cd", "cz", "cy", "cx", "cs", "cr", "cq", "cp", "cw", "cv", "cu", "ct", "pr", "ps", "pp", "pq", "pv", "pw", "pt", "pu", "pz", "px", "py", "wz", "pb", "pc", "pa", "pf", "pg", "pd", "pe", "pj", "pk", "ph", "pi", "pn", "po", "pl", "pm", "iy", "ix", "vb", "iz", "vd", "ve", "vf", "vg", "iq", "ip", "is", "ir", "iu", "it", "iw", "iv", "ii", "ih", "ik", "ij", "im", "il", "io", "in", "ia", "vy", "ic", "ib", "ie", "id", "ig", "if", "wh", "yi", "vr", "vs", "bd", "be", "bf", "bg", "ba", "bb", "bc", "bl", "bm", "bn", "bo", "bh", "bi", "bj", "bk", "bt", "bu", "bv", "bw", "bp", "bq", "br", "bs", "bx", "by", "bz", "oo", "on", "om", "ol", "ok", "oj", "oi", "oh", "og", "of", "oe", "od", "oc", "ob", "oa", "oz", "oy", "ox", "ow", "ov", "ou", "ot", "os", "or", "oq", "op", "hz", "hx", "hy", "hr", "hs", "hp", "hq", "hv", "hw", "ht", "hu", "hj", "hk", "hh", "hi", "hn", "ho", "hl", "hm", "hb", "hc", "ha", "hf", "hg", "hd", "he", "uy", "ux", "uz", "uu", "ut", "uw", "uv", "uq", "up", "us", "ur", "um", "ul", "uo", "un", "ui", "uh", "uk", "uj", "ue", "ud", "ug", "uf", "ua", "uc", "ub", "aa", "ac", "ab", "ae", "ad", "ag", "af", "ai", "ah", "ak", "aj", "am", "al", "ao", "an", "aq", "ap", "as", "ar", "au", "at", "aw", "av", "ay", "ax", "az", "nh", "ni", "nj", "nk", "nl", "nm", "nn", "no", "na", "nb", "nc", "nd", "ne", "nf", "ng", "nx", "ny", "nz", "np", "nq", "nr", "ns", "nt", "nu", "nv", "nw", "vp", "vq"};
    private static final int[] trans_vals = {0, 0, 3, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 4, 0, 15, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 11, 5, 0, 0, 0, 0, 3, 4, 0, 0, 0, 21, 49, 0, 0, 0, 4, 0, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 20, 0, 0, 0, 9, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 0, 0, 1, 0, 0, 0, 0, 16, 0, 0, 0, 0, 22, 0, 0, 3, 1, 0, 0, 15, 5, 1, 0, 2, 5, 0, 0, 4, 0, 0, 9, 0, 0, 1, 12, 20, 0, 0, 11, 0, 1, 0, 0, 0, 0, 0, 4, 9, 0, 1, 3, 0, 1, 3, 0, 0, 0, 7, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 2, 10, 1, 0, 0, 0, 0, 6, 21, 11, 16, 60, 0, 0, 0, 0, 5, 0, 0, 1, 4, 0, 2, 0, 0, 85, 0, 0, 0, 0, 2, 5, 29, 2, 10, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 24, 0, 3, 12, 0, 0, 2, 0, 7, 30, 0, 14, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 85, 1, 13, 0, 0, 15, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 7, 0, 5, 3, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 17, 0, 0, 0, 4, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 1, 0, 42, 13, 0, 35, 0, 6, 0, 0, 0, 0, 0, 9, 11, 5, 15, 0, 31, 8, 66, 3, 3, 1, 4, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 5, 0, 1, 0, 0, 5, 0, 0, 0, 4, 0, 2, 0, 5, 0, 0, 1, 7, 0, 0, 1, 1, 11, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 0, 0, 0, 15, 0, 0, 0, 0, 2, 0, 0, 0, 2, 11, 11, 1, 2, 20, 0, 2, 0, 0, 0, 1, 1, 2, 0, 22, 5, 0, 0, 2, 0, 1, 1, 0, 0, 19, 0, 1, 0, 4, 14, 10, 25, 0, 3, 3, 27, 31, 5, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 15, 0, 6, 2, 12, 0, 8, 0, 0, 0, 0, 0, 0, 6, 4, 0, 0, 0, 0, 0};

    private TreeMap<String, Integer> delMat = new TreeMap<>();
    private TreeMap<String, Integer> insMat = new TreeMap<>();
    private TreeMap<String, Integer> subMat = new TreeMap<>();
    private TreeMap<String, Integer> transMat = new TreeMap<>();

    private TreeMap<String, Integer> occurrences = new TreeMap<>();

    private int alphabetSize = 0;

    /**
     * Constructor
     * @param terms
     */
    public LetterProbability(ArrayList<String> terms) {
        buildMatrices();
        countOccurrences(terms);
    }

    /**
     * Build the actual confusion matrices out of the hard-coded arrays
     */
    private void buildMatrices() {
        for (int i = 0; i < del_keys.length; ++i) {
            delMat.put(del_keys[i], del_vals[i]);
        }

        for (int i = 0; i < ins_keys.length; ++i) {
            insMat.put(ins_keys[i], ins_vals[i]);
        }

        for (int i = 0; i < sub_keys.length; ++i) {
            subMat.put(sub_keys[i], sub_vals[i]);
        }

        for (int i = 0; i < trans_keys.length; ++i) {
            transMat.put(trans_keys[i], trans_vals[i]);
        }
    }

    /**
     * Count letter occurrences in the given list of terms, and the alphabet size.
     * @param terms
     */
    private void countOccurrences(ArrayList<String> terms) {
        for (String term : terms) {
            for (int i = 0; i < term.length() - 1; ++i) {
                String singleLetter = term.substring(i, i + 1);
                String doubleLetter = term.substring(i, i + 2);
                int singleCount = occurrences.getOrDefault(singleLetter, 0);
                int doubleCount = occurrences.getOrDefault(doubleLetter, 0);
                if (singleCount == 0) {  // Found a new single letter
                    ++alphabetSize;
                }
                occurrences.put(singleLetter, singleCount + 1);
                occurrences.put(doubleLetter, doubleCount + 1);
            }

            // Now take care of the last letter;
            String lastLetter = term.substring(term.length() - 1);
            int lastCount = occurrences.getOrDefault(lastLetter, 0);
            occurrences.put(lastLetter, lastCount + 1);
        }
    }

    //---------------------------------------------------//
    //********************* Getters *********************//
    //---------------------------------------------------//
    public TreeMap<String, Integer> getDelMat() {
        return delMat;
    }

    public TreeMap<String, Integer> getInsMat() {
        return insMat;
    }

    public TreeMap<String, Integer> getSubMat() {
        return subMat;
    }

    public TreeMap<String, Integer> getTransMat() {
        return transMat;
    }

    public TreeMap<String, Integer> getOccurrences() {
        return occurrences;
    }

    public int getAlphabetSize() {
        return alphabetSize;
    }
}
