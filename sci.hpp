//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
////                                                      ////
////    ███████╗ ██████╗██╗   ██╗  ██╗██████╗ ██████╗     ////
////    ██╔════╝██╔════╝██║   ██║  ██║██╔══██╗██╔══██╗    ////
////    ███████╗██║     ██║   ███████║██████╔╝██████╔╝    ////
////    ╚════██║██║     ██║   ██╔══██║██╔═══╝ ██╔═══╝     ////
////    ███████║╚██████╗██║██╗██║  ██║██║     ██║         ////
////    ╚══════╝ ╚═════╝╚═╝╚═╝╚═╝  ╚═╝╚═╝     ╚═╝         ////
////                                                      ////
////                (A cri5me™ production)                ////
////                                                      ////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



// TODO: Add BitSet type
// TODO: Maybe less dependence on standard 
//       headers? would be nice..
// TODO: Improve/Expand Allocators
//		  - realloc? it's quite unclear to me if realloc is really
//          worth it over freeing and allocating anew; maybe?
//          seems better to just allocate smarter lol (i.e. less frequently)
//        - alignment
// TODO: Add more string functions
//		  - strview -- I no longer recall what this was meant to do lol
//        - Actually, maybe upgrade str and friends...?
//           - Not sure what exactly I mean here lol; C++-ify slightly? *shrugs*
// TODO: Add (more) file I/O
//        - Actually a full virtual filesystem API would be nice
// TODO: Add fallback allocation with logging 
//		 (optional) to Temporary_Storage
// TODO: Add a simple logging API
//        - Just the basics: printf but to a "Log" which
//          can output to any combination of the console
//          and files on disk. Maybe we have other
//          options for output as well (i.e. a common interface)
//        - File rotation, ofc
// TODO: Add tmark and treset to Temporary_Storage?
// TODO: Add a debug memory allocator that can be used
// 		 as a wrapper around any Allocator* that provides
//		 some behind-the-scenes bookkeeping in order
//		 to detect things like leaks, double-frees, etc.
// TODO: Add an option for Arena*s to use guard pages
// TODO: Add an NBT-esque general-purpose binary data format.
//        - I don't want to just copy NBT, first of all. More betterer.
//          Would be nice if we could make API-side improvements...
//          I don't know. I wanted to use tinyrefl to use attributes to
//          generate serialization code from annotated structs, but tinyrefl
//          is a piece of shit so. And I haven't felt like writing a lexer and
//          making that work myself.. yet.
// TODO: Add string builder
//		  - maybe this uses the binary read/write thing?
//		  - or maybe the binary read/write thing uses this?
// TODO: Add hexdump printer
// TODO: Add a base64 codec
// TODO: Add RLE codec
// TODO: Add intel hex encoder
// 	 	  - Also SRec?
// TODO: Add sorting functions for Array<T>, etc.
// TODO: Add a segment tree
// TODO: Add a "Sparse Set" (the funky one w/ uninit memory)
//		  - Also, "Sparse Map", because it's free to add ^
// TODO: Add a ring buffer type
//        - Support both Static_Ring_Buffer and Ring_Buffer? or...
// TODO: Make sure this stuff works on winderps..... :/
// TODO: Add bignum structs
//		  - BigInt
//		  - BigRat?
//		  - BigDecimal?
//		  - Maybe also special things like u128/s128?
//			- May be a thing that already exists for larger
//			  CPU registers/SIMD *shrugs*
// TODO: Add a JSON parser/printer and "object model", as it were.
// TODO: Add intrinsics wrappers
// TODO: Add SIMD wrappers
//        - Probably largely based on structs with overloaded operators...
//          i.e. u32x4 or w/e
// TODO: Add various linked list structs
//        - intrusive
//        - non-intrusive? (is this even useful?)
//        - include iterators, etc.
// TODO: Add write/read varint to DataOutput/DataInput
// TODO: Add an LRU cache type
// TODO: Add a priority queue type
// TODO: Add command-line argument parsing



////////////////////////////////////
////////////////////////////////////
///            HEADER            ///
////////////////////////////////////
////////////////////////////////////



#ifndef SCI_H
#define SCI_H


/////////////////////
///    General    ///
/////////////////////


#include <stdlib.h>
#include <stdio.h> 
#include <stdarg.h>
#include <stddef.h>
#include <string.h>
#include <setjmp.h>
#include <assert.h>

#ifdef SCI_TESTING
#include <unistd.h>
#include <sys/wait.h>
#endif


#ifndef SCI_DEF
#define SCI_DEF extern
#endif


#ifndef __cplusplus
#error "Must compile as C++"
#endif


#if defined(_MSC_VER)
#define SCI_CC_MSVC
#elif defined(__clang__)
#define SCI_CC_CLANG
#elif defined(__GNUC__)
#define SCI_CC_GCC
#elif defined(__EMSCRIPTEN__)
#error "Emscripten not supported!"
#elif defined(__MINGW32__) || defined(__MINGW64__)
//#error "MinGW not yet supported!"
#else
#error "Unknown compiler!"
#endif


#if defined(__linux__)
#define SCI_OS_LINUX
#elif defined(_WIN32) || defined(_WIN64)
#define SCI_OS_WINDOWS
#elif defined(__APPLE__) && defined(__MACH__)
#define SCI_OS_OSX
#else
#warning "Unknown operating system!"
#endif


#define INTEGRAL_TYPES(X)      \
    X(u8, unsigned char)       \
    X(u16, unsigned short)     \
    X(u32, unsigned int)       \
    X(u64, unsigned long long) \
    X(s8, signed char)         \
    X(s16, signed short)       \
    X(s32, signed int)         \
    X(s64, signed long long)

#define FLOAT_TYPES(X) \
    X(f32, float)      \
    X(f64, double) 

#define X(name, ctype) using name = ctype;
INTEGRAL_TYPES(X)
FLOAT_TYPES(X)
#undef X


// Don't @ me.
static_assert(sizeof(u8)  == 1);
static_assert(sizeof(s8)  == 1);
static_assert(sizeof(u16) == 2);
static_assert(sizeof(s16) == 2);
static_assert(sizeof(u32) == 4);
static_assert(sizeof(s32) == 4);
static_assert(sizeof(u64) == 8);
static_assert(sizeof(s64) == 8);
static_assert(sizeof(f32) == 4);
static_assert(sizeof(f64) == 8);


#define CAT2(a, b) a##b
#define CAT(a, b) CAT2(a, b)


#define array_length(a) ((sizeof(a))/(sizeof(a[0])))
#define cast(t, v) ((t)(v))
#define unused(x) ((void)x)


#ifndef NULL
#define NULL 0
#endif


#ifndef offsetof
#define offsetof(type, member) ((u64)&(((type*)0)->member))
#endif


// TODO
#define nvrreturn 	  __attribute__((noreturn))
#define nvrinline 	  __attribute__((noinline))
#define forceinline	  __attribute__((always_inline))
#define static_init   __attribute__((constructor))
#define static_deinit __attribute__((destructor))


// TODO
#define debug_trap __builtin_trap
#define unreachable __builtin_unreachable


#ifndef PATH_SEPARATOR
#ifdef _WIN32
#define PATH_SEPARATOR '\\'
#else
#define PATH_SEPARATOR '/'
#endif
#endif


template<typename T>
void swap(T *a, T *b) {
	T c = *a;
	*a = *b;
	*b = c;
}

template<typename T>
void swap(T& a, T& b) {
	T c = a;
	a = b;
	b = c;
}


//////////////////////
///    Sections    ///
//////////////////////


#ifdef SCI_SECTIONS
#ifndef SCI_CC_GCC
#error "SCI_SECTIONS only supported on GCC"
#endif


#define SECTION(section_name) __attribute__((section(#section_name)))

#define SECTION_START(section_name) __start_##section_name[]
#define SECTION_STOP(section_name) __stop_##section_name[]

#define SECTION_START_SYMBOL(section_name, type)                            \
    ({                                                                      \
        extern const type SECTION_START(section_name);                      \
        __start_##section_name;                                             \
    })

#define SECTION_STOP_SYMBOL(section_name, type)                             \
    ({                                                                      \
        extern const type SECTION_STOP(section_name);                       \
        __stop_##section_name;                                              \
    })
    

#define SECTION_FOREACH(section_name, type, iter)                           \
    for(type const* iter = SECTION_START_SYMBOL(section_name, type);        \
        iter < SECTION_STOP_SYMBOL(section_name, type); iter++)
#endif


///////////////////
///    Maths    ///
///////////////////


// NOTE: I did this after trying to write each of these
// as a single function template taking type `t`. But yknow,
// sometimes C++ just doesn't feel like working so meh.
#define DEFROT(t)                             \
    constexpr t rotl(t x, t s) noexcept {     \
        t mask = sizeof(t) * 8 - 1;           \
        s &= mask;                            \
        return (x << s) | (x >> (-s & mask)); \
    }                                         \
    constexpr t rotr(t x, t s) noexcept {     \
        t mask = sizeof(t) * 8 - 1;           \
        s &= mask;                            \
        return (x >> s) | (x << (-s & mask)); \
    }
DEFROT(u8)
DEFROT(u16)
DEFROT(u32)
DEFROT(u64)
#undef DEFROT


// NOTE: And despite the NOTE above about the rotl/rotr functions,
// so far, is_pow2 has caused to trouble. WTF C++? Either work or don't!
template<typename T>
constexpr bool is_pow2(T x) noexcept { return x > 0 && (x & (x - 1)) == 0; }

constexpr u8 next_pow2(u8 y) noexcept {
	u8 x = y;
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x++;
	return x;
}

constexpr u16 next_pow2(u16 y) noexcept {
	u16 x = y;
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x++;
	return x;
}

constexpr u32 next_pow2(u32 y) noexcept {
	u32 x = y;
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	x++;
	return x;
}

constexpr u64 next_pow2(u64 y) noexcept {
	u64 x = y;
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	x |= x >> 32;
	x++;
	return x;
}

template<typename T>
constexpr T align_up(T x, T align) noexcept {
	assert(is_pow2(x) && "must align to a power of two");
    return (x + (align - 1)) & ~(align - 1);
}

template<typename T>
constexpr T align_down(T x, T align) noexcept {
	assert(is_pow2(x) && "must align to a power of two");
    return x - (x & (align - 1));
}

template<typename t>
constexpr t pow(t base, t exp) noexcept {
    if(base == 0 && exp == 0) {
        return 1;
    } else if(base == 0) {
        return 0;
    } else if(exp == 0) {
        return 1;
    }
    t result = 1;
    while(exp) {
        if(exp & 1) {
            result *= base;
        }
        exp >>= 1;
        base *= base;
    }
    return result;
}


#define BYTE_UNIT_MULTIPLES(X)  \
	X(1,  KILO,   KIBI) 		\
	X(2,  MEGA,   MEBI) 		\
	X(3,  GIGA,   GIBI) 		\
	X(4,  TERA,   TEBI) 		\
	X(5,  PETA,   PEBI) 		\
	X(6,  EXA,    EXI )

#define X(i, dec, bin) 						     \
	constexpr u64 dec##BYTE = pow<u64>(1000, i); \
	constexpr u64 bin##BYTE = pow<u64>(1024, i); \
	constexpr u64 dec##BYTES(u64 n) {		     \
		return n * dec##BYTE;				     \
	}										     \
	constexpr u64 bin##BYTES(u64 n) {		     \
		return n * bin##BYTE;				     \
	}
BYTE_UNIT_MULTIPLES(X)
#undef X


template<typename T>
constexpr T square(T x) noexcept { return x * x; }

template<typename T>
constexpr T sign(T x) noexcept { return x < 0 ? -1 : 1; }

template<typename T>
constexpr T clamp(T x, T min, T max) noexcept {
	if(x < min) return min;
	if(x > max) return max;
	return x;
}

template<typename T>
constexpr T min(T a, T b) noexcept {
	if(a < b) return a;
	return b;
}

template<typename T>
constexpr T max(T a, T b) noexcept {
	if(a > b) return a;
	return b;
}

template<typename T>
constexpr T abs(T x) noexcept {
	if(x < 0) return -x;
	return x;
}

template<typename T>
constexpr T lerp(T a, T b, T t) noexcept {
	return (1 - t) * a + t * b;
}

template<typename T>
constexpr T unlerp(T min, T max, T value) noexcept {
    return (value - min) / (max - min);
}

template<typename T>
constexpr T relerp(T in_min, T in_max, T value, T out_min, T out_max) noexcept {
    return lerp(out_min, out_max, unlerp(in_min, in_max, value));
}

template<typename T>
constexpr T move_towards(T value, T target, T rate) noexcept {
	if(abs(value - target) < 0.5) return target;
	return lerp(value, target, rate);
}


//////////////////////
///    OS Misc.    ///
//////////////////////


SCI_DEF u64 get_page_size();


///////////////////////
///    Allocator    ///
//////////////////////


struct Allocator {
	Allocator *parent;

    Allocator(Allocator *_parent = NULL) : parent(_parent) {}
    virtual void deinit() {}
    virtual void* alloc(u64 n) = 0;
    virtual void free(void *p) = 0;
};

SCI_DEF Allocator *allocator;

SCI_DEF void push_allocator(Allocator *a);
SCI_DEF void pop_allocator();
SCI_DEF void* xalloc(u64 n, Allocator *a = NULL);
SCI_DEF void xfree(void *p, Allocator *a = NULL);


struct _sci_new_wrapper{};
inline void* operator new(size_t, _sci_new_wrapper, void* ptr) { return ptr; }
inline void operator delete(void*, _sci_new_wrapper, void*) {}
inline void* operator new[](size_t, _sci_new_wrapper, void *ptr) { return ptr; }
inline void operator delete[](void*, _sci_new_wrapper, void*) {}
#define pnew(t, p, ...) (new(_sci_new_wrapper(), p) t(__VA_ARGS__))
#define xanew(t, a, ...) pnew(t, xalloc(sizeof(t), a), __VA_ARGS__)
#define xnew(t, ...) xanew(t, allocator, __VA_ARGS__)


///////////////////////////////
///    Temporary Storage    ///
///////////////////////////////


#ifndef TEMPORARY_STORAGE_SIZE
#define TEMPORARY_STORAGE_SIZE KIBIBYTES(64)
#endif

struct Temporary_Storage : public Allocator {
	u64 used;
    u64 high_water_mark;
    u8 data[TEMPORARY_STORAGE_SIZE];

    void* alloc(u64 n) override {
        assert(used + n < TEMPORARY_STORAGE_SIZE);
        void *p = &data[used];
        used += n;
        if(used > high_water_mark) high_water_mark = used;
        return p;
    }

    void free(void *p) override {}
};

SCI_DEF Allocator *temp_allocator;
SCI_DEF Temporary_Storage temporary_storage;

SCI_DEF void* talloc(u64 n);
SCI_DEF void treset();


///////////////////
///    Arena    ///
///////////////////


#ifndef ARENA_DEFAULT_BLOCK_SIZE
#define ARENA_DEFAULT_BLOCK_SIZE KIBIBYTES(64)
#endif

#ifndef ARENA_DEFAULT_PREALLOC_BLOCKS
#define ARENA_DEFAULT_PREALLOC_BLOCKS 1
#endif

#ifndef ARENA_MINIMUM_BLOCK_SIZE
#define ARENA_MINIMUM_BLOCK_SIZE KIBIBYTE
#endif

struct Arena : public Allocator {
	struct Block {
		struct Block *prev;
		u64 used = 0;
		u8 base[];

		Block(struct Block *_prev) : prev(_prev) {}
	};

	Block *current_block = NULL;
    Block *free_list = NULL;
	u64 block_size;

	Arena(Allocator *parent = allocator, u64 _block_size = ARENA_DEFAULT_BLOCK_SIZE, u32 prealloc_blocks = ARENA_DEFAULT_PREALLOC_BLOCKS) : Allocator(parent), block_size(_block_size) {
        for(u32 i = 0; i < prealloc_blocks; i++) {
            free_list = pnew(Block, xalloc(sizeof(Block) + block_size, parent), free_list);
        }
    }

    void deinit() override {
        auto b = current_block;
        while(b) {
            auto t = b->prev;
            xfree(b, parent);
            b = t;
        }
        b = free_list;
        while(b) {
            auto t = b->prev;
            xfree(b, parent);
            b = t;
        }
        Allocator::deinit();
    }

    void* alloc(u64 n) override {
        auto blk = current_block;

        if(!blk || blk->used + n >= block_size) {
            if(free_list) {
                blk = free_list;
                free_list = blk->prev;
            } else {
                blk = pnew(Block, xalloc(sizeof(Block) + block_size, parent), NULL);
            }
            current_block = blk;
        }

        assert(blk->used + n < block_size);

        auto p = blk->base + blk->used;
        blk->used += n;
        return p;
    }

    void free(void *p) override {}

    void clear() {
        auto b = current_block;
        while(b) {
            auto n = b->prev;
            b->prev = free_list;
            b->used = 0;
            free_list = b;
            b = n;
        }
        current_block = NULL;
    }
};


/////////////////////
///    Hashing    ///
/////////////////////


SCI_DEF u32 murmur3(void const *input, s32 len, u32 seed);
SCI_DEF u32 fnv1a(void const* input, u64 len);


/////////////////////
///    Strings    ///
/////////////////////


//
// NOTE: This is how our _str type works:
//
//    +-----------------+--------+------+
//    |  header (_str)  |  data  |  \0  |
//    +-----------------+--------+------+
//						|
//						+-> returned pointer
//

using str  = char*;       // str
using istr = char*;       // interned string
using cstr = char*;       // char*
using rstr = char const*; // char const*

struct _str {
	Allocator *a;
	u64 size;
	char data[];
};

#define strhdr(s) cast(_str*,(s)-sizeof(_str))

// TODO: use a hash table or something lmao
// TODO: allocate the _istrs linearly
struct _istr {
	_istr *prev;
	_str val;
};

#define istrhdr(s) cast(_istr*,(s)-sizeof(_istr))


SCI_DEF str mkstr(cstr s, u64 n, Allocator *a = NULL);
SCI_DEF str mkstr(cstr s, Allocator *a = NULL);
SCI_DEF void freestr(str s);
SCI_DEF u64 strsz(str s);
SCI_DEF str strcopy(str s, Allocator *a = NULL);
SCI_DEF bool streq(str a, str b);
SCI_DEF str substr(str s, u64 b, u64 e);
SCI_DEF str tvsprintf(rstr fmt, va_list args);
SCI_DEF str tsprintf(rstr fmt, ...);
SCI_DEF void tfprintf(FILE *fh, rstr fmt, ...);

SCI_DEF istr intern(cstr s);
SCI_DEF bool isintern(str s);


//////////////////////////
///    More General    ///
//////////////////////////


SCI_DEF void panic(rstr fmt, ...);
SCI_DEF void todo();


//////////////////////
///    File I/O    ///
//////////////////////


SCI_DEF str read_entire_file(str path, Allocator *a = NULL);
SCI_DEF bool write_entire_file(str path, str data);


//////////////////////////
///    Static Array    ///
//////////////////////////


template<typename T, u64 capacity>
struct Static_Array {
    static constexpr u64 size = capacity;
    using value_type = T;
	using type = Static_Array<T, capacity>;
    using ref_type = type&;
    using ref_const_type = type const&;
    using ptr_type = type*;
    using ptr_const_type = type const*;

    template<typename AT>
    struct IteratorBase {
        using type = IteratorBase<AT>;

        IteratorBase() : a(NULL) {}
        IteratorBase(AT _a, u64 _index) : a(_a), index(_index) {}

        type& operator++() { index++; return *this; }
        type operator++(int) { auto it = *this; operator++(); return it; }
        type& operator--() { index--; return *this; }
        type operator--(int) { auto it = *this; operator--(); return it; }

        T& operator*() { return (*a)[index]; }

        bool operator==(type const& b) const {
            if(a != b.a) return false;
            return index == b.index;
        }

        bool operator!=(type const& b) const { return !(*this == b ); }
    private:
        AT a;
        u64 index = 0;
    };

    using Iterator = IteratorBase<ptr_type>;
    using Const_Iterator = IteratorBase<ptr_const_type>;

	T data[capacity];
	u64 count = 0;

    Iterator begin() { return Iterator (this, 0); }
    Iterator end() { return Iterator(this, count); }
    Const_Iterator begin() const { return Const_Iterator(this, 0); }
    Const_Iterator end() const { return Const_Iterator(this, count); }

	void clear() {
		count = 0;
	}

	u64 push(T datum) {
		assert(count < capacity);
		data[count] = datum;
		return count++;
	}

	T pop() {
		assert(count > 0);
		return data[--count];
	}

	T& operator[](u64 x) { return data[x]; }
	T const& operator[](u64 x) const { return data[x]; }
};


////////////////////////////
///    Slot Allocator    ///
////////////////////////////


// NOTE: This is extremely similar to Static_Array
// but, at the moment, I feel like keeping them separate
// since they have different intended uses.
//          - sci4me, 5/12/20

template<typename T, u32 size>
struct Slot_Allocator {
    T slots[size];
    s32 count;

    Slot_Allocator() {
        memset(slots, 0, sizeof(slots));
        count = 0;
    }

    void clear() {
        count = 0;
    }

    s32 index_of(T x) {
        for(s32 i = 0; i < count; i++) {
            if(slots[i] == x) return i;
        }
        return -1;
    }

    s32 alloc(T x) {      
        s32 i = index_of(x);
        if(i != -1) return i;

        if(count == size) return -1;

        slots[count] = x;
        return count++;
    }
};


///////////////////
///    Array    ///
///////////////////


#ifndef ARRAY_DEFAULT_SIZE
#define ARRAY_DEFAULT_SIZE 64
#endif

template<typename T>
struct Array {
    struct Iterator {
        Iterator() : a(NULL) {}
        Iterator(Array<T> *_a, u64 _index) : a(_a), index(_index) {}

        Iterator& operator++() { index++; return *this; }
        Iterator operator++(int) { auto it = *this; operator++(); return it; }
        Iterator& operator--() { index--; return *this; }
        Iterator operator--(int) { auto it = *this; operator--(); return it; }

        T& operator*() { return (*a)[index]; }

        bool operator==(Iterator const& b) const {
            if(a != b.a) return false;
            return index == b.index;
        }

        bool operator!=(Iterator const& b) const { return !(*this == b ); }
    private:
        Array<T> *a;
        u64 index = 0;
    };

    struct Const_Iterator {
        Const_Iterator() : a(NULL) {}
        Const_Iterator(Array<T> const* _a, u64 _index) : a(_a), index(_index) {}

        Const_Iterator& operator++() { index++; return *this; }
        Const_Iterator operator++(int) { auto it = *this; operator++(); return it; }
        Const_Iterator& operator--() { index--; return *this; }
        Const_Iterator operator--(int) { auto it = *this; operator--(); return it; }

        T const& operator*() { return (*a)[index]; }

        bool operator==(Const_Iterator const& b) const {
            if(a != b.a) return false;
            return index == b.index;
        }

        bool operator!=(Const_Iterator const& b) const { return !(*this == b ); }
    private:
        Array<T> const* a;
        u64 index = 0;
    };

	u64 count = 0;
	u64 size = 0;
	T *data = NULL;
	Allocator *a;

	Array(Allocator *_a = ::allocator) : a(_a) {}

    Iterator begin() { return Iterator(this, 0); }
    Iterator end() { return Iterator(this, count); }
    Const_Iterator begin() const { return Const_Iterator(this, 0); }
    Const_Iterator end() const { return Const_Iterator(this, count); }

	void free() {
		if(data) {
			xfree(data, a);
			count = 0;
			size = 0;
			data = NULL;
		}
	}

	void resize(u64 size) {
		if(size < ARRAY_DEFAULT_SIZE) return;

		assert(size >= count);

		T *new_data = cast(T*, xalloc(sizeof(T) * size, a));
		assert(new_data);
		memcpy(new_data, data, sizeof(T) * count);
		xfree(data, a);

		this->size = size;
		this->data = new_data;
	}

	void clear() {
		count = 0;
	}

	u64 push(T v) {
		check_init();

		if(count == size) {
			resize(size * 2);
		}

		u64 i = count++;
		data[i] = v;
		return i;
	}

	T pop() {
		assert(data);
		assert(count > 0);
		return data[--count];
	}

	void insert(T v, u64 index) {
		check_init();

		if(count == size) {
			resize(size * 2);
		}

		for(u64 i = count - 1; i >= index; i--) {
			if(i > count) break;
			data[i + 1] = data[i];
		}
		data[index] = v;
		count++;
	}

	T ordered_remove(u64 i) {
		assert(data);
		assert(i < count);
		T r = data[i];
		for(u64 j = i; j < count - 1; j++) {
			data[j] = data[j + 1];
		}
		count--;
		return r;
	}

	T unordered_remove(u64 i) {
		assert(data);
		assert(i < count);
		T r = data[i];
		data[i] = data[count - 1];
		count--;
		return r;
	}

	bool contains(T x) {
		if(!data) return false;
		return index(x) != -1;
	}

	s64 index(T x) {
		for(u64 i = 0; i < count; i++) {
			if(data[i] == x) return i;
		}
		return -1;
	}

	T& operator[](u64 x) { return data[x]; }
	T const& operator[](u64 x) const { return data[x]; }

private:
	void check_init() {
		if(data == NULL) {
			size = ARRAY_DEFAULT_SIZE;
			data = cast(T*, xalloc(sizeof(T) * size, a));
		}
	}
};


////////////////////////
///    Hash Table    ///
////////////////////////


// TODO: use 64-bit hashes?
// TODO: switch hash table to use a 64-bit size and count?


// NOTE: We're currently just using a fixed seed
// theoretically we _could_ generate it randomly
// at app-startup (or we could even be more granular
// than that and make it unique to the hash table, but
// meh, don't know that we need to.)
// I just got this number off of random.org.
//				- sci4me, 5/21/20
constexpr u32 HASH_TABLE_DEFAULT_SEED = 0xB23D66D5;

template<typename K>
using HashFN = u32 (*)(K const&); // TODO: Better name?

template<typename K>
using EqFN = bool (*)(K const&, K const&); // TODO: Better name?

template<typename T>
u32 default_hash_fn(T const& v) {
	return murmur3((void const *) &v, sizeof(T), HASH_TABLE_DEFAULT_SEED);
}

template<typename T>
bool default_eq_fn(T const& a, T const& b) {
	return memcmp((void const*) &a, (void const*) &b, sizeof(T)) == 0;
}


#ifndef HASH_TABLE_LOAD_FACTOR_SHRINK_THRESHOLD
#define HASH_TABLE_LOAD_FACTOR_SHRINK_THRESHOLD 0.1f
#endif
#ifndef HASH_TABLE_LOAD_FACTOR_EXPAND_THRESHOLD
#define HASH_TABLE_LOAD_FACTOR_EXPAND_THRESHOLD 0.85f
#endif
#ifndef HASH_TABLE_DEFAULT_SIZE
#define HASH_TABLE_DEFAULT_SIZE 256
#endif


// NOTE: Currently, we _require_ size to be a power of 2!
// Eventually, we _should_ switch to using sizes that are
// prime numbers. Probably.
// 				- sci4me, 5/21/20

// NOTE: This hash table uses "Robin Hood Hashing":
// " Robin Hood hashing is a technique for implementing hash tables.   				"
// " It is based on open addressing with a simple but clever twist: As new   		"
// " keys are inserted, old keys are shifted around in a way such that all   		"
// " keys stay reasonably close to the slot they originally hash to. In particular, "
// " the variance of the keys distances from their "home" slots is minimized.  		"
//				- sci4me, 11/23/20

template<typename K, typename V, HashFN<K> hash_fn = default_hash_fn, EqFN<K> eq_fn = default_eq_fn>
struct Hash_Table {
	struct Slot {
		K key;
		V value;
		u32 hash;
	};

	u32 count;
	u32 size;
	u32 mask;
	Slot *slots;

	// TODO: make this a lazy-init structure
	void init(u32 size = HASH_TABLE_DEFAULT_SIZE) {
		assert(is_pow2(size));
		this->count = 0;
		this->size = size;
		this->mask = size - 1;
		slots = cast(Slot*, xalloc(size * sizeof(Slot)));
		for(u32 i = 0; i < size; i++) slots[i].hash = 0;
	}

	void free() {
		xfree(slots);
	}

	void resize(u32 new_size) {
		assert(new_size);
		new_size = next_pow2(new_size);
        assert(is_pow2(new_size));
        assert(new_size >= count);

		u32 old_count = count;
		u32 old_size = size;
		Slot *old_slots = slots;

		count = 0;
		size = new_size;
		mask = size - 1;
		slots = cast(Slot*, xalloc(size * sizeof(Slot)));
		for(u32 i = 0; i < size; i++) slots[i].hash = 0;

		for(u32 i = 0; i < old_size; i++) {
			if(old_slots[i].hash) {
				set(old_slots[i].key, old_slots[i].value);
			}
		}

		assert(count == old_count);

		xfree(old_slots);
	}

	void set(K _key, V _value) {
		if(load_factor() > HASH_TABLE_LOAD_FACTOR_EXPAND_THRESHOLD) {
			resize(size * 2);
		}

		K key = _key;
		V value = _value;
		u32 hash = hash_key(key);

		s32 i = hash & mask;
		s32 dist = 0;
		for(;;) {
			if(slots[i].hash == 0) {
				slots[i].key = key;
				slots[i].value = value;
				slots[i].hash = hash;
				count++;
				return;
			}

            s32 epd = probe_dist(i, slots[i].hash);
			if(epd < dist) {
				assert(slots[i].hash);

				K _k = slots[i].key;
				V _v = slots[i].value;
				u32 _h = slots[i].hash;

				slots[i].key = key;
				slots[i].value = value;
				slots[i].hash = hash;

				key = _k;
				value = _v;
				hash = _h;

				dist = epd;
			}

			i = (i + 1) & mask;
			dist++;
		}
	}

	V get(K key) const {
		s32 i = index_of(key);
		if(i == -1) {
			V dummy;
			memset(&dummy, 0, sizeof(V)); // NOTE: not strictly necessary...
			return dummy;
		}
		return slots[i].value;
	}

	s32 index_of(K key) const {
		u32 hash = hash_key(key);
		s32 i = hash & mask;
		u32 dist = 0;
		for(;;) {
			if(slots[i].hash == 0) {
				return -1;
			}

            // NOTE TODO: Are we supposed to be using `hash`?
            // Or are we supposed to be using slots[i].hash?...
            s32 epd = probe_dist(i, hash);
			if(dist > epd) {
				return -1;
			}

			if(slots[i].hash == hash && eq_fn(slots[i].key, key)) {
				return i;
			}

			i = (i + 1) & mask;
			dist++;
		}
	}

	bool remove(K key) {
		s32 i = index_of(key);
		if(i == -1) return false;

		for(s32 j = 0; j < size - 1; j++) {
			s32 k = (i + 1) & mask;

			if(slots[k].hash == 0) break;

            s32 epd = probe_dist(k, slots[k].hash);
			if(epd == 0) break;

			memcpy(&slots[i], &slots[k], sizeof(Slot));

			i = k;
		}

		slots[i].hash = 0;
		count--;

		if(load_factor() < HASH_TABLE_LOAD_FACTOR_SHRINK_THRESHOLD) {
			resize(max(size / 2, u32(HASH_TABLE_DEFAULT_SIZE)));
		}

		return true;
	}

	f32 load_factor() const {
		return (f32) count / (f32) size;
	}

private:
    s32 probe_dist(u32 hash, u32 slot_index) const noexcept {
        return (slot_index + size - (hash & mask)) & mask;
    }

	static u32 hash_key(K key) {
		u32 h = hash_fn(key);
		// NOTE: a hash of 0 represents an empty slot
		if(h == 0) h |= 1;
		return h;
	}
};


u32 str_hash_fn(str const& s) {
	return murmur3((void const*) s, strsz(s), HASH_TABLE_DEFAULT_SEED);
}

bool str_eq_fn(str const& a, str const& b) {
	return streq(a, b);
}

u32 cstr_hash_fn(cstr const& s) {
    return murmur3((void const*) s, strlen(s), HASH_TABLE_DEFAULT_SEED);
}

bool cstr_eq_fn(cstr const& a, cstr const& b) {
    return strcmp(a, b) == 0;
}


/////////////////////////
///    I/O Streams    ///
/////////////////////////


struct DataInput {
    virtual u8 read_u8() = 0;

    inline u16 read_u16() {
        return (cast(u16, read_u8()) << 8) | 
                cast(u16, read_u8());
    }

    inline u32 read_u32() {
        return (cast(u32, read_u8()) << 24) |
               (cast(u32, read_u8()) << 16) |
               (cast(u32, read_u8()) << 8) |
                read_u8();
    }

    inline u64 read_u64() {
        return (cast(u64, read_u8()) << 56) |
               (cast(u64, read_u8()) << 48) |
               (cast(u64, read_u8()) << 40) |
               (cast(u64, read_u8()) << 32) |
               (cast(u64, read_u8()) << 24) |
               (cast(u64, read_u8()) << 16) |
               (cast(u64, read_u8()) << 8) |
                cast(u64, read_u8());
    }

    inline s8 read_s8() { return cast(s8, read_u8()); }
    inline s16 read_s16() { return cast(s16, read_u16()); }
    inline s32 read_s32() { return cast(s32, read_u32()); }
    inline s64 read_s64() { return cast(s64, read_u64()); }

    inline f32 read_f32() {
        union {
            u32 u;
            f32 f;
        };
        u = read_u32();
        return f;
    }

    inline f64 read_f64() {
        union {
            u64 u;
            f64 f;
        };
        u = read_u64();
        return f;
    }

    inline str read_str() {
        u64 len = read_u64();
        str s = mkstr(NULL, len, allocator);
        for(u64 i = 0; i < len; i++) s[i] = read_u8();
        return s;
    }

    inline cstr read_cstr() {
        u64 len = read_u64();
        cstr s = cast(cstr, xalloc(len + 1, allocator));
        for(u64 i = 0; i < len; i++) s[i] = read_u8();
        return s;
    }
};


struct DataOutput {
    virtual void write_u8(u8 x) = 0;

    inline void write_u16(u16 x) {
        write_u8((x >> 8) & 0xFF);
        write_u8(x & 0xFF);
    }

    inline  void write_u32(u32 x) {
        write_u8((x >> 24) & 0xFF);
        write_u8((x >> 16) & 0xFF);
        write_u8((x >> 8) & 0xFF);
        write_u8(x & 0xFF);
    }

    inline void write_u64(u64 x) {
        write_u8((x >> 56) & 0xFF);
        write_u8((x >> 48) & 0xFF);
        write_u8((x >> 40) & 0xFF);
        write_u8((x >> 32) & 0xFF);
        write_u8((x >> 24) & 0xFF);
        write_u8((x >> 16) & 0xFF);
        write_u8((x >> 8) & 0xFF);
        write_u8(x & 0xFF);
    }

    inline void write_s8(s8 x) { write_u8(cast(u8, x)); }
    inline void write_s16(s16 x) { write_u16(cast(u16, x)); }
    inline void write_s32(s32 x) { write_u32(cast(u32, x)); }
    inline void write_s64(s64 x) { write_u64(cast(u64, x)); }

    inline void write_f32(f32 x) {
        union {
            f32 f;
            u32 u;
        };
        f = x;
        write_u32(u);
    }

    inline void write_f64(f64 x) {
        union {
            f64 f;
            u64 u;
        };
        f = x;
        write_u64(u); 
    }

    inline void write_str(str s) {
        u64 len = strsz(s);
        write_u64(len);
        for(u64 i = 0; i < len; i++) write_u8(s[i]);
    }

    inline void write_cstr(cstr s) {
        u64 len = strlen(s);
        write_u64(len);
        for(u64 i = 0; i < len; i++) write_u8(s[i]);
    }
};

struct RandomAccessDataOutput {
    virtual u64 reserve(u64 n) = 0;
    virtual void set_u8(u64 i, u8 x) = 0;

    inline void set_u16(u64 i, u16 x) {
        set_u8(i, (x >> 8) & 0xFF);
        set_u8(i + 1, x & 0xFF);
    }

    inline void set_u32(u64 i, u32 x) {
        set_u8(i, (x >> 24) & 0xFF);
        set_u8(i + 1, (x >> 16) & 0xFF);
        set_u8(i + 2, (x >> 8) & 0xFF);
        set_u8(i + 3, x & 0xFF);
    }

    inline void set_u64(u64 i, u64 x) {
        set_u8(i, (x >> 56) & 0xFF);
        set_u8(i + 1, (x >> 48) & 0xFF);
        set_u8(i + 2, (x >> 40) & 0xFF);
        set_u8(i + 3, (x >> 32) & 0xFF);
        set_u8(i + 4, (x >> 24) & 0xFF);
        set_u8(i + 5, (x >> 16) & 0xFF);
        set_u8(i + 6, (x >> 8) & 0xFF);
        set_u8(i + 7, x & 0xFF);
    }

    inline void set_s8(u64 i, s8 x) { set_u8(i, cast(u8, x)); }
    inline void set_s16(u64 i, s16 x) { set_u16(i, cast(u16, x)); }
    inline void set_s32(u64 i, s32 x) { set_u32(i, cast(u32, x)); }
    inline void set_s64(u64 i, s64 x) { set_u64(i, cast(u64, x)); }
};


#ifndef BYTEBUF_DEFAULT_SIZE
#define BYTEBUF_DEFAULT_SIZE 4096
#endif

struct ByteBuf : public DataInput, public DataOutput, public RandomAccessDataOutput {
    static ByteBuf wrap(cstr s) { return ByteBuf(cast(u8*, s), 0, strlen(s)); }
    static ByteBuf wrap(u8 *p, u64 n) { return ByteBuf(p, 0, n); }

    u8 *data = NULL;
    u64 index = 0;
    u64 size;

    ByteBuf() : ByteBuf(BYTEBUF_DEFAULT_SIZE) {}
    ByteBuf(u64 _size) : size(_size) {}
    ByteBuf(u8 *_data, u64 _index, u64 _size) : data(_data), index(_index), size(_size) {}

    void deinit() {
        xfree(data);
    }

    void write_to_file(cstr path) {
        FILE *fh = fopen(path, "w");
        assert(fwrite(data, sizeof(char), index, fh) == index);
        fflush(fh);
        fclose(fh);
    }

    bool read_from_file(cstr path) {
        FILE *fh = fopen(path, "r");
        if(!fh) return false;

        fseek(fh, 0L, SEEK_END);
	    u64 size = ftell(fh);
	    rewind(fh);

        data = cast(u8*, xalloc(size));
    	assert(fread(data, sizeof(char), size, fh) == size);

        fclose(fh);

        return true;
    }

    str tostr() {
        auto s = mkstr(NULL, index);
        memcpy(s, data, index);
        return s;
    }

    void reset() {
        index = 0;
    }

    void expand() {
        if(data) {
            u64 new_size = size * 2;
            u8 *new_data = cast(u8*, xalloc(new_size));
            memcpy(new_data, data, size);
            xfree(data);
            data = new_data;
            size = new_size;
        } else {
            assert(size);
            data = cast(u8*, xalloc(size));
        }
    }

    u8 read_u8() override {
        return data[index++];
    }

    void write_u8(u8 x) override {
        if(index >= size || data == NULL) {
            expand();
        }
    
        data[index++] = x;
    }

    u64 reserve(u64 n) override {
        if(index + n >= size || data == NULL) {
            expand();
        }
        u64 i = index;
        index += n;
        return i;
    }

    void set_u8(u64 i, u8 x) override {
        data[i] = x;
    }
};


/////////////////////
///    Testing    ///
/////////////////////


#ifdef SCI_TESTING

#ifndef SCI_SECTIONS
// TODO: Probably don't require this since it depends on a GCC extension...
#error "SCI_TESTING requires SCI_SECTIONS"
#endif

struct Test {
    void (*func)();
    rstr name;
	rstr file;
	u8 pad[4];
};

#define deftest(func_name)                              \
    void test_##func_name();                            \
    const Test SECTION(v2cc_tests)  	 				\
    test_info_##func_name = {                           \
        .func = test_##func_name,                       \
        .name = #func_name,                             \
		.file = __FILE__								\
    };                                                  \
    void test_##func_name()

#define TESTS_FOREACH(iter) SECTION_FOREACH(v2cc_tests, Test, iter)

#else

#define deftest(func_name) void test_##func_name()

#define TESTS_FOREACH(iter) if(0)

#endif


SCI_DEF s32 run_tests();


#undef SCI_DEF

#endif



////////////////////////////////////
////////////////////////////////////
///        IMPLEMENTATION        ///
////////////////////////////////////
////////////////////////////////////



#if defined(SCI_IMPL) && !defined(SCI_IMPL_DONE)
#define SCI_IMPL_DONE 


#define SCI_MANGLE(x) _sci_hpp_internal_##x


//////////////////////
///    OS Misc.    ///
//////////////////////


#ifdef SCI_OS_WINDOWS
#define WIN32_LEAN_AND_MEAN
#include "windows.h"
#endif


static u64 SCI_MANGLE(page_size);

u64 get_page_size() { return SCI_MANGLE(page_size); }

static_init void init_page_size() {
#ifdef SCI_OS_WINDOWS
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    SCI_MANGLE(page_size) = cast(u64, info.dwPageSize);
#else
    SCI_MANGLE(page_size) = sysconf(_SC_PAGESIZE);
#endif
}


///////////////////////
///    Allocator    ///
///////////////////////


struct Sys_Allocator : public Allocator {
    void* alloc(u64 n) override {
        void *p = malloc(n);
        memset(p, 0, n);
        return p;
    }

    void free(void *p) override {
        ::free(p);
    }
};


Sys_Allocator sys_allocator;
Allocator *allocator = &sys_allocator;


void push_allocator(Allocator *a) {
    todo();
}

void pop_allocator() {
    todo();
}

void* xalloc(u64 n, Allocator *a) {
	auto a2 = a ? a : allocator;
	return a2->alloc(n);
}

void xfree(void *p, Allocator *a) {
	auto a2 = a ? a : allocator;
	a2->free(p);
}


///////////////////////////////
///    Temporary Storage    ///
///////////////////////////////


Allocator *temp_allocator;
Temporary_Storage temporary_storage;

static_init void temp_alloc_init() {
    memset(&temporary_storage, 0, sizeof(Temporary_Storage));
    temp_allocator = &temporary_storage;
}


void* talloc(u64 n) {
	return temporary_storage.alloc(n);
}

void treset() {
	temporary_storage.used = 0;
}


/////////////////////
///    Hashing    ///
/////////////////////


u32 murmur3(void const *input, s32 len, u32 seed) {
	constexpr u32 C1 = 0xCC9E2D51;
	constexpr u32 C2 = 0x1B873593;

	u8 const *data = (u8 const*) input;
	s32 n_blocks = len / 4;

	u32 h = seed;

	u32 const *blocks = (u32 const*) (data + n_blocks * 4);
	for(s32 i = -n_blocks; i; i++) {
		u32 k = blocks[i];

		k *= C1;
		k = rotl(k, 15);
		k *= C2;

		h ^= k;
		h = rotl(h, 13);
		h = h * 5 + 0xE6546B64;
	}

	u8 const *tail = (u8 const*) (data + n_blocks * 4);
	u32 k = 0;
	switch(len & 3) {
		case 3: k ^= tail[2] << 16; [[fallthrough]];
		case 2: k ^= tail[1] << 8;  [[fallthrough]];
		case 1: k ^= tail[0];
				k *= C1;
				k = rotl(k, 15);
				k *= C2;
				h ^= k;
	}

	h ^= len;

	h ^= h >> 16;
	h *= 0x85EBCA6B;
	h ^= h >> 13;
	h *= 0xC2B2AE35;
	h ^= h >> 16;

	return h;
}


u32 fnv1a(void const* input, u64 len) {
	u32 x = 0x811C9DC5;
	for(u64 i = 0; i < len; i++) {
		x ^= cast(u8*, input)[i];
		x *= 0x01000193;
	}
	return x;
}


/////////////////////
///    Strings    ///
/////////////////////


str mkstr(cstr s, u64 n, Allocator *a) {
	auto a2 = a ? a : allocator;
	_str *r = cast(_str*, xalloc(sizeof(_str) + n + 1, a2));
	r->a = a2;
	r->size = n;
	if(s) memcpy(r->data, s, n);
	r->data[n] = 0;
	return r->data;
}

str mkstr(cstr s, Allocator *a) {
	return mkstr(s, strlen(cast(rstr, s)), a);
}

void freestr(str s) {
	auto s2 = strhdr(s);
	s2->a->free(s2);
}

u64 strsz(str s) {
	return strhdr(s)->size;
}

str strcopy(str s, Allocator *a) {
	auto s2 = strhdr(s);
	auto a2 = a ? a : s2->a;
	_str *r = cast(_str*, a2->alloc(sizeof(_str) + s2->size + 1));
	r->a = a2;
	r->size = s2->size;
	memcpy(r->data, s2->data, r->size + 1);
	return r->data;
}

bool streq(str a, str b) {
	if(strsz(a) != strsz(b)) return false;
	return strcmp(a, b) == 0;
}

str substr(str s, u64 b, u64 e) {
	auto s2 = strhdr(s);
	assert(e >= b);
	assert(e < s2->size);
	u64 n = e - b + 1;
	return mkstr(s2->data + b, n);
}

str tvsprintf(rstr fmt, va_list args) {
	va_list args2;
	va_copy(args2, args);
	u64 n = vsnprintf(NULL, 0, fmt, args2);
	va_end(args2);
	auto r = cast(_str*, xalloc(sizeof(_str) + n + 1, temp_allocator));
	vsprintf(r->data, fmt, args);
	r->data[n] = 0;
    r->size = n;
	return r->data;
}

str tsprintf(rstr fmt, ...) {
	va_list args;
	va_start(args, fmt);
	return tvsprintf(fmt, args);
}

void tfprintf(FILE *fh, rstr fmt, ...) {
	va_list args;
	va_start(args, fmt);
	str s = tvsprintf(fmt, args);
	va_end(args);
	fprintf(fh, "%s", s);
}


static _istr *interned_strings = NULL;

static_deinit void free_interned_strings() {
	auto cur = interned_strings;
	while(cur) {
		auto prev = cur->prev;
		xfree(cur);
		cur = prev;
	}
}

istr intern(cstr s) {
	// NOTE: used rather than strsz to support regular C strings
	u64 n = strlen(cast(rstr, s));

	auto cur = interned_strings;
	while(cur) {
		auto is = cast(istr, cur->val.data);
		if(n == strlen(is) && memcmp(s, is, n) == 0) {
			return is;
		}
		cur = cur->prev;
	}

	auto r = cast(_istr*, xalloc(sizeof(_istr) + n + 1));
	r->prev = interned_strings;
	r->val.a = allocator;
	r->val.size = n;
	memcpy(r->val.data, s, n);
	r->val.data[n] = 0;
	interned_strings = r;
	return r->val.data;
}

bool isintern(str s) {
	auto cur = interned_strings;
	while(cur) {
		auto is = cast(istr, cur->val.data);
		if(s == is) return true;
		cur = cur->prev;
	}
	return false;
}


//////////////////////////
///    More General    ///
//////////////////////////


nvrreturn void panic(rstr fmt, ...) {
	va_list args;
	va_start(args, fmt);
	str s = tvsprintf(fmt, args);
	va_end(args);
	tfprintf(stderr, "panic: %s\n\n", s);
	*((volatile u32*)0) = 42;
	exit(EXIT_FAILURE);
}


nvrreturn void todo() {
	panic("TODO");
}


//////////////////////
///    File I/O    ///
//////////////////////


str read_entire_file(str path, Allocator *a) {
	FILE *fh = fopen(path, "rb");
	if(!fh) return NULL;

	fseek(fh, 0L, SEEK_END);
	u64 size = ftell(fh);
	rewind(fh);

	_str *r = cast(_str*, xalloc(sizeof(_str) + size + 1, a));
	r->a = allocator;
	r->size = size;

	assert(fread(r->data, sizeof(char), size, fh) == size);
	r->data[size] = 0;

	fclose(fh);

	return cast(str, r->data);
}

bool write_entire_file(str path, str data) {
	FILE *fh = fopen(path, "wb");
	if(!fh) return false;

	u64 n = strsz(data);
	assert(fwrite(data, sizeof(char), n, fh) == n);

	fflush(fh);
	fclose(fh);

	return true;
}


/////////////////////
///    Testing    ///
/////////////////////


s32 run_tests() {
#ifdef SCI_TESTING
	u32 passed = 0;
    u32 failed = 0;
    u32 test_count = 0;

    TESTS_FOREACH(t) {
        test_count++;

        auto pid = fork();
        if(pid == 0) {
            t->func();
            _exit(0);
        } else {
            int status;
            pid_t cpid;
            assert((cpid = wait(&status)) == pid);
            if(status) {
                printf(" \u001b[31;1m*\u001b[0m %s:%s\n", basename(t->file), t->name);
                failed++;
            } else {
                printf(" \u001b[32;1m\u2713\u001b[0m %s:%s\n", basename(t->file), t->name);
                passed++;
            }
        }
    }    

    printf(
		"\nFailed: %u (%0.2f%%)\nPassed: %u (%0.2f%%)\nTotal:  %u\n", 
		failed, 
		((f64)failed / (f64)test_count) * 100.0f, 
		passed, 
		((f64)passed / (f64)test_count) * 100.0f, 
		test_count
	);

    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
#else
	fprintf(stderr, "compiled without SCI_TESTING\n");
	return EXIT_FAILURE;
#endif
}



#undef SCI_MANGLE


#endif



// MIT License
// 
// Copyright (c) 2021 Scitoshi Nakayobro
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.