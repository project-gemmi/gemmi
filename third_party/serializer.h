#ifndef ZPP_SERIALIZER_H
#define ZPP_SERIALIZER_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#if __cplusplus >= 201703L
#include <optional>
#include <variant>
#endif
#ifdef ZPP_SERIALIZER_FREESTANDING
#include <string_view>
#else
#include <mutex>
#include <shared_mutex>
#include <unordered_map>
#endif

namespace zpp
{
/**
 * Supports serialization of objects and polymorphic objects.
 * Example of non polymorphic serialization:
 * ~~~
 * class point
 * {
 * public:
 *     point() = default;
 *     point(int x, int y) noexcept :
 *           m_x(x),
 *           m_y(y)
 *       {
 *       }
 *
 *     friend zpp::serializer::access;
 *     template <typename Archive, typename Self>
 *     static void serialize(Archive & archive, Self & self)
 *     {
 *         archive(self.m_x, self.m_y);
 *     }
 *
 *     int get_x() const noexcept
 *     {
 *           return m_x;
 *       }
 *
 *       int get_y() const noexcept
 *       {
 *           return m_y;
 *       }
 *
 * private:
 *     int m_x = 0;
 *     int m_y = 0;
 * };
 *
 * static void foo()
 * {
 *     std::vector<unsigned char> data;
 *     zpp::serializer::memory_input_archive in(data);
 *     zpp::serializer::memory_output_archive out(data);
 *
 *     out(point(1337, 1338));
 *
 *     point my_point;
 *     in(my_point);
 *
 *     std::cout << my_point.get_x() << ' ' << my_point.get_y() << '\n';
 * }
 * ~~~
 *
 * Example of polymorphic serialization:
 * ~~~
 * class person : public zpp::serializer::polymorphic
 * {
 * public:
 *     person() = default;
 *     explicit person(std::string name) noexcept :
 *           m_name(std::move(name))
 *       {
 *       }
 *
 *     friend zpp::serializer::access;
 *     template <typename Archive, typename Self>
 *     static void serialize(Archive & archive, Self & self)
 *     {
 *           archive(self.m_name);
 *       }
 *
 *       const std::string & get_name() const noexcept
 *       {
 *           return m_name;
 *       }
 *
 *     virtual void print() const
 *     {
 *         std::cout << "person: " << m_name;
 *     }
 *
 * private:
 *     std::string m_name;
 * };
 *
 * class student : public person
 * {
 * public:
 *     student() = default;
 *     student(std::string name, std::string university) noexcept :
 *         person(std::move(name)),
 *         m_university(std::move(university))
 *     {
 *     }
 *
 *     friend zpp::serializer::access;
 *     template <typename Archive, typename Self>
 *     static void serialize(Archive & archive, Self & self)
 *     {
 *         person::serialize(archive, self);
 *           archive(self.m_university);
 *       }
 *
 *     virtual void print() const
 *     {
 *         std::cout << "student: " << person::get_name() << ' ' <<
 * m_university << '\n';
 *     }
 *
 * private:
 *     std::string m_university;
 * };
 *
 * namespace
 * {
 * zpp::serializer::register_types<
 *    zpp::serializer::make_type<person,
 * zpp::serializer::make_id("v1::person")>,
 *    zpp::serializer::make_type<student,
 * zpp::serializer::make_id("v1::student")> > _; } // <anynymous namespace>
 *
 * static void foo()
 * {
 *     std::vector<unsigned char> data;
 *     zpp::serializer::memory_input_archive in(data);
 *     zpp::serializer::memory_output_archive out(data);
 *
 *     std::unique_ptr<person> my_person =
 * std::make_unique<student>("1337", "1337University"); out(my_person);
 *
 *     my_person = nullptr;
 *     in(my_person);
 *
 *     my_person->print();
 * }
 *
 * static void bar()
 * {
 *     std::vector<unsigned char> data;
 *     zpp::serializer::memory_input_archive in(data);
 *     zpp::serializer::memory_output_archive out(data);
 *
 *     out(zpp::serializer::as_polymorphic(student("1337",
 * "1337University")));
 *
 *     std::unique_ptr<person> my_person;
 *     in(my_person);
 *
 *     my_person->print();
 * }
 * ~~~
 */
namespace serializer
{
#ifdef ZPP_SERIALIZER_FREESTANDING
namespace freestanding
{
/**
 * Returns the error category for a given error code enumeration type,
 * using an argument dependent lookup of a user implemented category
 * function.
 */
template <typename ErrorCode>
decltype(auto) category()
{
    return category(ErrorCode{});
}

/**
 * The error category which responsible for translating error codes to
 * error messages.
 */
class error_category
{
public:
    /**
     * Returns the error category name.
     */
    virtual std::string_view name() const noexcept = 0;

    /**
     * Return the error message for a given error code.
     * For success codes, it is unspecified what value is returned.
     * For convienience, you may return zpp::error::no_error for success.
     * All other codes must return non empty string views.
     */
    virtual std::string_view message(int code) const noexcept = 0;

    /**
     * Returns true if success code, else false.
     */
    bool success(int code) const
    {
        return code == m_success_code;
    }

protected:
    /**
     * Creates an error category whose success code is 'success_code'.
     */
    constexpr error_category(int success_code) :
        m_success_code(success_code)
    {
    }

    /**
     * Destroys the error category.
     */
    ~error_category() = default;

private:
    /**
     * The success code.
     */
    int m_success_code{};
};

/**
 * Creates an error category, whose name and success
 * code are specified, as well as a message translation
 * logic that returns the error message for every error code.
 * Note: message translation must not throw.
 */
template <typename ErrorCode, typename Messages>
constexpr auto make_error_category(std::string_view name,
                                   ErrorCode success_code,
                                   Messages && messages)
{
    // Create a category with the name and messages.
    class category final : public error_category,
                           private std::remove_reference_t<Messages>
    {
    public:
        constexpr category(std::string_view name,
                           ErrorCode success_code,
                           Messages && messages) :
            error_category(
                std::underlying_type_t<ErrorCode>(success_code)),
            std::remove_reference_t<Messages>(
                std::forward<Messages>(messages)),
            m_name(name)
        {
        }

        std::string_view name() const noexcept override
        {
            return m_name;
        }

        std::string_view message(int code) const noexcept override
        {
            return this->operator()(ErrorCode{code});
        }

    private:
        std::string_view m_name;
    } category(name, success_code, std::forward<Messages>(messages));

    // Return the category.
    return category;
}

/**
 * Represents an error to be initialized from an error code
 * enumeration.
 * The error code enumeration must have 'int' as underlying type.
 * Defining an error code enum and a category for it goes as follows.
 * Example:
 * ~~~
 * namespace my_namespace
 * {
 * enum class my_error : int
 * {
 *     success = 0,
 *     something_bad = 1,
 *     something_really_bad = 2,
 * };
 *
 * inline const zpp::error_category & category(my_error)
 * {
 *     constexpr static auto error_category =
 *         zpp::make_error_category("my_category", my_error::success,
 *             [](auto code) -> std::string_view {
 *                 switch (code) {
 *                     case my_error::success:
 *                         return zpp::error::no_error;
 *                     case my_code:something_bad:
 *                         return "Something bad happened.";
 *                     case my_error::something_really_bad:
 *                         return "Something really bad happened.";
 *                     default:
 *                         return "Unknown error occurred.";
 *                 }
 *             }
 *         );
 *     return error_category;
 * }
 * } // my_namespace
 * ~~~
 */
class error
{
public:
    /**
     * Disables default construction.
     */
    error() = delete;

    /**
     * Constructs an error from an error code enumeration, the
     * category is looked by using argument dependent lookup on a
     * function named 'category' that receives the error code
     * enumeration value.
     */
    template <typename ErrorCode>
    error(ErrorCode error_code) :
        m_category(std::addressof(
            zpp::serializer::freestanding::category<ErrorCode>())),
        m_code(std::underlying_type_t<ErrorCode>(error_code))
    {
    }

    /**
     * Constructs an error from an error code enumeration, the
     * category is given explicitly in this overload.
     */
    template <typename ErrorCode>
    error(ErrorCode error_code, const error_category & category) :
        m_category(std::addressof(category)),
        m_code(std::underlying_type_t<ErrorCode>(error_code))
    {
    }

    /**
     * Returns the error category.
     */
    const error_category & category() const
    {
        return *m_category;
    }

    /**
     * Returns the error code.
     */
    int code() const
    {
        return m_code;
    }

    /**
     * Returns the error message. Calling this on
     * a success error is implementation defined according
     * to the error category.
     */
    std::string_view message() const
    {
        return m_category->message(m_code);
    }

    /**
     * Returns true if the error indicates success, else false.
     */
    explicit operator bool() const
    {
        return m_category->success(m_code);
    }

    /**
     * No error message value.
     */
    static constexpr std::string_view no_error{};

private:
    /**
     * The error category.
     */
    const error_category * m_category{};

    /**
     * The error code.
     */
    int m_code{};
};
} // namespace freestanding

enum class error : int
{
    success = 0,
    out_of_range = 1,
    variant_is_valueless = 2,
    null_pointer_serialization = 3,
};

inline const freestanding::error_category & category(error)
{
    constexpr static auto error_category =
        freestanding::make_error_category(
            "zpp::serializer",
            error::success,
            [](auto code) -> std::string_view {
                switch (code) {
                case error::success:
                    return freestanding::error::no_error;
                case error::out_of_range:
                    return "[zpp::serializer] Out of range error";
                case error::variant_is_valueless:
                    return "[zpp::serializer] Cannot serialize a "
                           "valueless variant.";
                case error::null_pointer_serialization:
                    return "[zpp::serializer] Cannot serialize a null "
                           "pointer.";
                default:
                    return "[zpp::serializer] Unknown error occurred.";
                }
            });
    return error_category;
}
#endif // ZPP_SERIALIZER_FREESTANDING

namespace detail
{
/**
 * Map any sequence of types to void.
 */
template <typename...>
using void_t = void;

/**
 * Tests if all conditions are true, empty means true.
 * Example:
 * ~~~
 * all_of<true, false, true>::value == false
 * all_of<true, true>::value == true
 * all_of<false, false>::value == false
 * all_of<>::value == true
 * ~~~
 */
template <bool... bConditions>
struct all_of : std::true_type
{
};

template <bool... bConditions>
struct all_of<false, bConditions...> : std::false_type
{
};

template <bool... bConditions>
struct all_of<true, bConditions...> : all_of<bConditions...>
{
};

template <>
struct all_of<true> : std::true_type
{
};

/**
 * Remove const of container value_type
 */
template <typename Container, typename = void>
struct container_nonconst_value_type
{
    using type = std::remove_const_t<typename Container::value_type>;
};

/**
 * Same as above, except in case of std::map and std::unordered_map, and
 * similar, we also need to remove the const of the key type.
 */
template <template <typename...> class Container,
          typename KeyType,
          typename MappedType,
          typename... ExtraTypes>
struct container_nonconst_value_type<
    Container<KeyType, MappedType, ExtraTypes...>,
    void_t<
        // Require existence of key_type.
        typename Container<KeyType, MappedType, ExtraTypes...>::key_type,

        // Require existence of mapped_type.
        typename Container<KeyType, MappedType, ExtraTypes...>::
            mapped_type,

        // Require that the value type is a pair of const KeyType and
        // MappedType.
        std::enable_if_t<std::is_same<
            std::pair<const KeyType, MappedType>,
            typename Container<KeyType, MappedType, ExtraTypes...>::
                value_type>::value>>>
{
    using type = std::pair<KeyType, MappedType>;
};

/**
 * Alias to the above.
 */
template <typename Container>
using container_nonconst_value_type_t =
    typename container_nonconst_value_type<Container>::type;

/**
 * The serializer exception template.
 */
template <typename Base, std::size_t Id>
class exception : public Base
{
public:
    /**
     * Use the constructors from the base class.
     */
    using Base::Base;
};

/**
 * A no operation, single byte has same representation in little/big
 * endian.
 */
inline constexpr std::uint8_t swap_byte_order(std::uint8_t value) noexcept
{
    return value;
}

/**
 * Swaps the byte order of a given integer.
 */
inline constexpr std::uint16_t
swap_byte_order(std::uint16_t value) noexcept
{
    return (std::uint16_t(swap_byte_order(std::uint8_t(value))) << 8) |
           (swap_byte_order(std::uint8_t(value >> 8)));
}

/**
 * Swaps the byte order of a given integer.
 */
inline constexpr std::uint32_t
swap_byte_order(std::uint32_t value) noexcept
{
    return (std::uint32_t(swap_byte_order(std::uint16_t(value))) << 16) |
           (swap_byte_order(std::uint16_t(value >> 16)));
}

/**
 * Swaps the byte order of a given integer.
 */
inline constexpr std::uint64_t
swap_byte_order(std::uint64_t value) noexcept
{
    return (std::uint64_t(swap_byte_order(std::uint32_t(value))) << 32) |
           (swap_byte_order(std::uint32_t(value >> 32)));
}

/**
 * Rotates the given number left by count bits.
 */
template <typename Integer>
constexpr auto rotate_left(Integer number, std::size_t count)
{
    return (number << count) | (number >> ((sizeof(number) * 8) - count));
}

/**
 * Checks if has 'data()' member function.
 */
template <typename Type, typename = void>
struct has_data_member_function : std::false_type
{
};

/**
 * Checks if has 'data()' member function.
 */
template <typename Type>
struct has_data_member_function<
    Type,
    void_t<decltype(std::declval<Type &>().data())>> : std::true_type
{
};

} // namespace detail

/**
 * @name Exceptions
 * @{
 */
using out_of_range = detail::exception<std::out_of_range, 0>;
using undeclared_polymorphic_type_error =
    detail::exception<std::runtime_error, 1>;
using attempt_to_serialize_null_pointer_error =
    detail::exception<std::logic_error, 2>;
using polymorphic_type_mismatch_error =
    detail::exception<std::runtime_error, 3>;
using attempt_to_serialize_valueless_variant =
    detail::exception<std::runtime_error, 4>;
using variant_index_out_of_range = detail::exception<std::out_of_range, 5>;
/**
 * @}
 */

/**
 * If C++17 or greater, use shared mutex, else, use shared timed mutex.
 */
#ifndef ZPP_SERIALIZER_FREESTANDING
#if __cplusplus >= 201703L
/**
 * The shared mutex type, defined to shared mutex when available.
 */
using shared_mutex = std::shared_mutex;
#else
/**
 * The shared mutex type, defined to shared timed mutex when shared mutex
 * is not available.
 */
using shared_mutex = std::shared_timed_mutex;
#endif
#endif

/**
 * The base class for polymorphic serialization.
 */
class polymorphic
{
public:
    /**
     * Pure virtual destructor, in order to become abstract
     * and make derived classes polymorphic.
     */
    virtual ~polymorphic() = 0;
};

/**
 * Default implementation for the destructor.
 */
inline polymorphic::~polymorphic() = default;

/**
 * Allow serialization with saving (output) archives, of objects held by
 * reference, that will be serialized as polymorphic, meaning, with leading
 * polymorphic serialization id.
 */
template <typename Type>
class polymorphic_wrapper
{
public:
    /**
     * Constructs from the given object to be serialized as polymorphic.
     */
    explicit polymorphic_wrapper(const Type & object) noexcept :
        m_object(object)
    {
        static_assert(std::is_base_of<polymorphic, Type>::value,
                      "The given type is not derived from polymorphic");
    }

    /**
     * Returns the object to be serialized as polymorphic.
     */
    const Type & operator*() const noexcept
    {
        return m_object;
    }

private:
    /**
     * The object to be serialized as polymorphic.
     */
    const Type & m_object;
}; // polymorphic_wrapper

/**
 * A facility to save object with leading polymorphic serialization id.
 */
template <typename Type>
auto as_polymorphic(const Type & object) noexcept
{
    return polymorphic_wrapper<Type>(object);
}

/**
 * The size type of the serializer.
 * It is used to indicate the size for containers.
 */
using size_type = std::uint32_t;

/**
 * The serialization id type,
 */
using id_type = std::uint64_t;

/**
 * This class grants the serializer access to the serialized types.
 */
class access
{
public:
    /**
     * Allows placement construction of types.
     */
    template <typename Item, typename... Arguments>
    static auto
    placement_new(void * pAddress, Arguments &&... arguments) noexcept(
        noexcept(Item(std::forward<Arguments>(arguments)...)))
    {
        return ::new (pAddress)
            Item(std::forward<Arguments>(arguments)...);
    }

    /**
     * Allows dynamic construction of types.
     * This overload is for non polymorphic serialization.
     */
    template <typename Item,
              typename... Arguments,
              typename = std::enable_if_t<
                  !std::is_base_of<polymorphic, Item>::value>>
    static auto make_unique(Arguments &&... arguments)
    {
        // Construct the requested type, using new since constructor might
        // be private.
        return std::unique_ptr<Item>(
            new Item(std::forward<Arguments>(arguments)...));
    }

    /**
     * Allows dynamic construction of types.
     * This overload is for polymorphic serialization.
     */
    template <typename Item,
              typename... Arguments,
              typename = std::enable_if_t<
                  std::is_base_of<polymorphic, Item>::value>,
              typename = void>
    static auto make_unique(Arguments &&... arguments)
    {
        // We create a deleter that will delete using the base class
        // polymorphic which, as we declared, has a public virtual
        // destructor.
        struct deleter
        {
            void operator()(Item * item) noexcept
            {
                delete static_cast<polymorphic *>(item);
            }
        };

        // Construct the requested type, using new since constructor might
        // be private.
        return std::unique_ptr<Item, deleter>(
            new Item(std::forward<Arguments>(arguments)...));
    }

    /**
     * Allows destruction of types.
     */
    template <typename Item>
    static void destruct(Item & item) noexcept
    {
        item.~Item();
    }
}; // access

/**
 * Enables serialization of arbitrary byte data.
 * Use only with care.
 */
template <typename Item>
class bytes
{
public:
    /**
     * Constructs the bytes wrapper from pointer and count of items.
     */
    bytes(Item * items, std::size_t count) : m_items(items), m_count(count)
    {
    }

    /**
     * Returns a pointer to the first item.
     */
    Item * data() const noexcept
    {
        return m_items;
    }

    /**
     * Returns the size in bytes of the bytes data.
     */
    std::size_t size_in_bytes() const noexcept
    {
        return m_count * sizeof(Item);
    }

    /**
     * Returns the count of items in the bytes wrapper.
     */
    std::size_t count() const noexcept
    {
        return m_count;
    }

private:
    /**
     * Pointer to the items.
     */
    Item * m_items = nullptr;

    /**
     * The number of items.
     */
    std::size_t m_count{};
};

/**
 * Allows serialization as bytes data.
 * Use only with care.
 */
template <typename Item>
bytes<Item> as_bytes(Item * item, std::size_t count)
{
    static_assert(std::is_trivially_copyable<Item>::value,
                  "Must be trivially copyable");

    return {item, count};
}

/**
 * Allows serialization as bytes data.
 * Use only with care.
 */
inline bytes<unsigned char> as_bytes(void * data, std::size_t size)
{
    return {static_cast<unsigned char *>(data), size};
}

/**
 * Allows serialization as bytes data.
 */
inline bytes<const unsigned char> as_bytes(const void * data,
                                           std::size_t size)
{
    return {static_cast<const unsigned char *>(data), size};
}

/**
 * The serialization method type.
 */
template <typename Archive, typename = void>
struct serialization_method;

/**
 * The serialization method type exporter, for loading (input) archives.
 */
template <typename Archive>
struct serialization_method<Archive, typename Archive::loading>
{
    /**
     * Disabled default constructor.
     */
    serialization_method() = delete;

    /**
     * The exported type.
     */
    using type = void (*)(Archive &, std::unique_ptr<polymorphic> &);
}; // serialization_method

/**
 * The serialization method type exporter, for saving (output) archives.
 */
template <typename Archive>
struct serialization_method<Archive, typename Archive::saving>
{
    /**
     * Disabled default constructor.
     */
    serialization_method() = delete;

    /**
     * The exported type.
     */
    using type = void (*)(Archive &, const polymorphic &);
}; // serialization_method

/**
 * The serialization method type.
 */
template <typename Archive>
using serialization_method_t =
    typename serialization_method<Archive>::type;

/**
 * Make a serialization method from type and a loading (input) archive.
 */
template <typename Archive,
          typename Type,
          typename...,
          typename = typename Archive::loading>
serialization_method_t<Archive> make_serialization_method() noexcept
{
    return [](Archive & archive, std::unique_ptr<polymorphic> & object) {
        auto concrete_type = access::make_unique<Type>();
        archive(*concrete_type);
        object.reset(concrete_type.release());
    };
}

/**
 * Make a serialization method from type and a saving (output) archive.
 */
template <typename Archive,
          typename Type,
          typename...,
          typename = typename Archive::saving,
          typename = void>
serialization_method_t<Archive> make_serialization_method() noexcept
{
    return [](Archive & archive, const polymorphic & object) {
        archive(dynamic_cast<const Type &>(object));
    };
}

/**
 * This is the base archive of the serializer.
 * It enables saving and loading items into/from the archive, via
 * operator().
 */
template <typename ArchiveType>
class archive
{
public:
    /**
     * The derived archive type.
     */
    using archive_type = ArchiveType;

    /**
     * Save/Load the given items into/from the archive.
     */
    template <typename... Items>
    auto operator()(Items &&... items)
    {
        // Disallow serialization of pointer types.
        static_assert(
            detail::all_of<!std::is_pointer<
                std::remove_reference_t<Items>>::value...>::value,
            "Serialization of pointer types is not allowed");

        // Serialize the items.
        return serialize_items(std::forward<Items>(items)...);
    }

protected:
    /**
     * Constructs the archive.
     */
    archive() = default;

    /**
     * Protected destructor to allow safe public inheritance.
     */
    ~archive() = default;

private:
    /**
     * Serialize the given items, one by one.
     */
    template <typename Item, typename... Items>
    auto serialize_items(Item && first, Items &&... items)
    {
#ifndef ZPP_SERIALIZER_FREESTANDING
        // Invoke serialize_item the first item.
        serialize_item(std::forward<Item>(first));
#else
        // Invoke serialize_item the first item.
        if (auto result = serialize_item(std::forward<Item>(first));
            !result) {
            return result;
        }
#endif
        // Serialize the rest of the items.
        return serialize_items(std::forward<Items>(items)...);
    }

    /**
     * Serializes zero items.
     */
    auto serialize_items()
    {
#ifdef ZPP_SERIALIZER_FREESTANDING
        return freestanding::error{error::success};
#endif
    }

    /**
     * Serialize a single item.
     * This overload is for class type items with serialize method.
     */
    template <typename Item,
              typename...,
              typename = decltype(std::remove_reference_t<Item>::serialize(
                  std::declval<archive_type &>(), std::declval<Item &>()))>
    auto serialize_item(Item && item)
    {
        // Forward as lvalue.
        return std::remove_reference_t<Item>::serialize(concrete_archive(),
                                                        item);
    }

    /**
     * Serialize a single item.
     * This overload is for types with outer serialize method.
     */
    template <typename Item,
              typename...,
              typename = decltype(serialize(std::declval<archive_type &>(),
                                            std::declval<Item &>())),
              typename = void>
    auto serialize_item(Item && item)
    {
        // Forward as lvalue.
        return serialize(concrete_archive(), item);
    }

    /**
     * Serialize a single item.
     * This overload is for fundamental types.
     */
    template <typename Item,
              typename...,
              typename = std::enable_if_t<std::is_fundamental<
                  std::remove_reference_t<Item>>::value>,
              typename = void,
              typename = void>
    auto serialize_item(Item && item)
    {
        // Forward as lvalue.
        return concrete_archive().serialize(item);
    }

    /**
     * Serialize a single item.
     * This overload is for enum classes.
     */
    template <typename Item,
              typename...,
              typename = std::enable_if_t<
                  std::is_enum<std::remove_reference_t<Item>>::value>,
              typename = void,
              typename = void,
              typename = void>
    auto serialize_item(Item && item)
    {
        // If the enum is const, we want the type to be a const type, else
        // non-const.
        using integral_type = std::conditional_t<
            std::is_const<std::remove_reference_t<Item>>::value,
            const std::underlying_type_t<std::remove_reference_t<Item>>,
            std::underlying_type_t<std::remove_reference_t<Item>>>;

        // Cast the enum to the underlying type, and forward as lvalue.
        return concrete_archive().serialize(
            reinterpret_cast<std::add_lvalue_reference_t<integral_type>>(
                item));
    }

    /**
     * Serialize bytes data.
     */
    template <typename T>
    auto serialize_item(bytes<T> && item)
    {
        return concrete_archive().serialize(item.data(),
                                            item.size_in_bytes());
    }

    /**
     * Returns the concrete archive.
     */
    archive_type & concrete_archive()
    {
        return static_cast<archive_type &>(*this);
    }
}; // archive

/**
 * This archive serves as an output archive, which saves data into memory.
 * Every save operation appends data into the vector or view type.
 * This archive serves as an optimization and type erasure for the view
 * type for polymorphic serialization.
 */
class basic_memory_output_archive
    : public archive<basic_memory_output_archive>
{
public:
    /**
     * The base archive.
     */
    using base = archive<basic_memory_output_archive>;

    /**
     * Declare base as friend.
     */
    friend base;

    /**
     * saving archive.
     */
    using saving = void;

protected:
    /**
     * Constructs a memory output archive, that outputs to the given
     * vector.
     */
    explicit basic_memory_output_archive(
        std::vector<unsigned char> & output) noexcept :
        m_output_vector(std::addressof(output))
    {
    }

    /**
     * Constructs a memory output archive, that outputs to the given
     * view.
     */
    basic_memory_output_archive(unsigned char * data,
                                std::size_t size) noexcept :
        m_data(data),
        m_capacity(size)
    {
    }

    /**
     * Serialize a single item - save its data.
     */
    template <typename Item>
    auto serialize(Item && item)
    {
        // Check if we are about to go beyond the capacity.
        if (m_offset + sizeof(item) > m_capacity) {
            if (!m_output_vector) {
#ifndef ZPP_SERIALIZER_FREESTANDING
                throw out_of_range(
                    "Serialization to view type archive is out of range.");
#else
                return freestanding::error{error::out_of_range};
#endif
            }
            m_capacity = (m_capacity + sizeof(item)) * 3 / 2;
            m_output_vector->resize(m_capacity);
            m_data = m_output_vector->data();
        }

        // Copy the data to the end of the view.
        std::copy_n(
            reinterpret_cast<const unsigned char *>(std::addressof(item)),
            sizeof(item),
            m_data + m_offset);

        // Increase the offset.
        m_offset += sizeof(item);

#ifdef ZPP_SERIALIZER_FREESTANDING
        return freestanding::error{error::success};
#endif
    }

    /**
     * Serialize bytes data - save its data.
     */
    auto serialize(const void * data, std::size_t size)
    {
        // Check if we are about to go beyond the capacity.
        if (m_offset + size > m_capacity) {
            if (!m_output_vector) {
#ifndef ZPP_SERIALIZER_FREESTANDING
                throw out_of_range(
                    "Serialization to view type archive is out of range.");
#else
                return freestanding::error{error::out_of_range};
#endif
            }
            m_capacity = (m_capacity + size) * 3 / 2;
            m_output_vector->resize(m_capacity);
            m_data = m_output_vector->data();
        }

        // Copy the data to the end of the output.
        std::copy_n(static_cast<const unsigned char *>(data),
                    size,
                    m_data + m_offset);

        // Increase the offset.
        m_offset += size;

#ifdef ZPP_SERIALIZER_FREESTANDING
        return freestanding::error{error::success};
#endif
    }

    /**
     * Resizes the vector to the desired size.
     */
    void fit_vector()
    {
        m_output_vector->resize(m_offset);
    }

    /**
     * Refresh the vector.
     */
    void refresh_vector() noexcept
    {
        m_data = m_output_vector->data();
        m_capacity = m_output_vector->size();
        m_offset = m_capacity;
    }

    /**
     * Returns the data pointer.
     */
    unsigned char * data() const noexcept
    {
        return m_data;
    }

    /**
     * Returns the current offset.
     */
    std::size_t offset() const noexcept
    {
        return m_offset;
    }

    /**
     * Returns the current offset.
     */
    void reset(std::size_t offset = {}) noexcept
    {
        m_offset = offset;
    }

private:
    /**
     * The output vector, may be null in which case working with
     * view type represented by data and size below.
     */
    std::vector<unsigned char> * m_output_vector{};

    /**
     * Points to the data.
     */
    unsigned char * m_data{};

    /**
     * The current capacity size of the data.
     */
    std::size_t m_capacity{};

    /**
     * The offset of the output data.
     */
    std::size_t m_offset{};
}; // basic_memory_output_archive

/**
 * This archive serves as an output archive, which saves data into memory.
 * Every save operation appends data into a user provided view type.
 */
class memory_view_output_archive : private basic_memory_output_archive
{
public:
    /**
     * The base archive.
     */
    using base = basic_memory_output_archive;

    /**
     * Constructing the view from pointer and size.
     */
    memory_view_output_archive(unsigned char * data,
                               std::size_t size) noexcept :
        basic_memory_output_archive(data, size)
    {
    }

    /**
     * Constructs a memory output archive, that outputs to the given
     * view type.
     */
    template <typename View,
              typename...,
              typename ViewType = std::remove_reference_t<View>,
              typename = decltype(std::declval<View &>().size()),
              typename = decltype(std::declval<View &>().begin()),
              typename = decltype(std::declval<View &>().end()),
              typename = decltype(std::declval<View &>().data()),
              typename = std::enable_if_t<
                  std::is_trivially_destructible<ViewType>::value>,
              typename = std::enable_if_t<
                  std::is_same<typename ViewType::value_type,
                               unsigned char>::value &&
                  std::is_base_of<std::random_access_iterator_tag,
                                  typename std::iterator_traits<
                                      typename ViewType::iterator>::
                                      iterator_category>::value>>
    explicit memory_view_output_archive(View && view) noexcept :
        memory_view_output_archive(view.data(), view.size())
    {
    }

    /**
     * Saves items into the archive.
     */
    template <typename... Items>
    auto operator()(Items &&... items)
    {
        // Save previous offset.
        auto offset = this->offset();

#ifndef ZPP_SERIALIZER_FREESTANDING
        try {
            // Serialize the items.
            base::operator()(std::forward<Items>(items)...);
        } catch (...) {
            this->reset(offset);
            throw;
        }
#else
        // Serialize the items.
        auto result = base::operator()(std::forward<Items>(items)...);
        if (!result) {
            this->reset(offset);
        }
        return result;
#endif
    }

    /**
     * Returns the data pointer.
     */
    using base::data;

    /**
     * Returns the current offset in the data.
     */
    using base::offset;

    /**
     * Allow to reset offset for advanced use.
     */
    using base::reset;
};

/**
 * This archive serves as an output archive, which saves data into memory.
 * Every save operation appends data into the vector.
 */
class memory_output_archive : private basic_memory_output_archive
{
public:
    /**
     * The base archive.
     */
    using base = basic_memory_output_archive;

    /**
     * Constructs a memory output archive, that outputs to the given
     * vector.
     */
    explicit memory_output_archive(
        std::vector<unsigned char> & output) noexcept :
        basic_memory_output_archive(output)
    {
    }

    /**
     * Saves items into the archive.
     */
    template <typename... Items>
    auto operator()(Items &&... items)
    {
        refresh_vector();

        // The original offset.
        auto offset = this->offset();

#ifndef ZPP_SERIALIZER_FREESTANDING
        // Serialize the items.
        try {
            base::operator()(std::forward<Items>(items)...);
            fit_vector();
        } catch (...) {
            this->reset(offset);
            throw;
        }
#else
        // Serialize the items.
        auto result = base::operator()(std::forward<Items>(items)...);
        if (!result) {
            this->reset(offset);
            return result;
        }

        fit_vector();
        return result;
#endif
    }

    /**
     * Returns the data pointer.
     */
    using base::data;

    /**
     * Returns the current offset in the data.
     */
    using base::offset;

    /**
     * Allow to reset offset for advanced use.
     */
    using base::reset;
};

/**
 * This archive serves as the memory view input archive, which loads data
 * from non owning memory. Every load operation advances an offset to that
 * the next data may be loaded on the next iteration.
 */
class memory_view_input_archive : public archive<memory_view_input_archive>
{
public:
    /**
     * The base archive.
     */
    using base = archive<memory_view_input_archive>;

    /**
     * Declare base as friend.
     */
    friend base;

    /**
     * Loading archive.
     */
    using loading = void;

    /**
     * Construct a memory view input archive, that loads data from an array
     * of given pointer and size.
     */
    memory_view_input_archive(const unsigned char * input,
                              std::size_t size) noexcept :
        m_input(input),
        m_size(size)
    {
    }

    /**
     * Constructs a memory view input  archive, that loads data from
     * the given view type.
     */
    template <typename View,
              typename...,
              typename ViewType = std::remove_reference_t<View>,
              typename = decltype(std::declval<View &>().size()),
              typename = decltype(std::declval<View &>().begin()),
              typename = decltype(std::declval<View &>().end()),
              typename = decltype(std::declval<View &>().data()),
              typename = std::enable_if_t<
                  std::is_trivially_destructible<ViewType>::value>,
              typename = std::enable_if_t<
                  std::is_same<
                      std::remove_const_t<typename ViewType::value_type>,
                      unsigned char>::value &&
                  std::is_base_of<std::random_access_iterator_tag,
                                  typename std::iterator_traits<
                                      typename ViewType::iterator>::
                                      iterator_category>::value>>
    explicit memory_view_input_archive(View && view) noexcept :
        memory_view_input_archive(view.data(), view.size())
    {
    }

    /**
     * Returns the current offset in the input data.
     */
    const unsigned char * data() const noexcept
    {
        return m_input;
    }

    /**
     * Returns the current offset in the input data.
     */
    std::size_t offset() const noexcept
    {
        return m_offset;
    }

    /**
     * Resets the serialization to offset, to allow advanced use.
     */
    void reset(std::size_t offset = {}) noexcept
    {
        m_offset = offset;
    }

protected:
    /**
     * Serialize a single item - load it from the vector.
     */
    template <typename Item>
    auto serialize(Item && item)
    {
        // Verify that the vector is large enough to contain the item.
        if (m_size < (sizeof(item) + m_offset)) {
#ifndef ZPP_SERIALIZER_FREESTANDING
            throw out_of_range("Input vector was not large enough to "
                               "contain the requested item");
#else
            return freestanding::error{error::out_of_range};
#endif
        }

        // Fetch the item from the vector.
        std::copy_n(
            m_input + m_offset,
            sizeof(item),
            reinterpret_cast<unsigned char *>(std::addressof(item)));

        // Increase the offset according to item size.
        m_offset += sizeof(item);

#ifdef ZPP_SERIALIZER_FREESTANDING
        return freestanding::error{error::success};
#endif
    }

    /**
     * Serializes bytes data.
     */
    auto serialize(void * data, std::size_t size)
    {
        // Verify that the vector is large enough to contain the data.
        if (m_size < (size + m_offset)) {
#ifndef ZPP_SERIALIZER_FREESTANDING
            throw out_of_range("Input vector was not large enough to "
                               "contain the requested item");
#else
            return freestanding::error{error::out_of_range};
#endif
        }

        // Fetch the bytes data from the vector.
        std::copy_n(
            m_input + m_offset, size, static_cast<unsigned char *>(data));

        // Increase the offset according to data size.
        m_offset += size;

#ifdef ZPP_SERIALIZER_FREESTANDING
        return freestanding::error{error::success};
#endif
    }

private:
    /**
     * The input data.
     */
    const unsigned char * m_input{};

    /**
     * The input size.
     */
    std::size_t m_size{};

    /**
     * The next input.
     */
    std::size_t m_offset{};
}; // memory_view_input_archive

/**
 * This archive serves as the memory input archive, which loads data from
 * owning memory. Every load operation erases data from the beginning of
 * the vector.
 */
class memory_input_archive : private memory_view_input_archive
{
public:
    /**
     * The base archive.
     */
    using base = memory_view_input_archive;

    /**
     * Construct a memory input archive from a vector.
     */
    memory_input_archive(std::vector<unsigned char> & input) :
        memory_view_input_archive(input.data(), input.size()),
        m_input(std::addressof(input))
    {
    }

    /**
     * Load items from the archive.
     */
    template <typename... Items>
    auto operator()(Items &&... items)
    {
        // Update the input archive.
        static_cast<memory_view_input_archive &>(*this) = {
            m_input->data(), m_input->size()};

        // Save the original offset.
        auto offset = this->offset();

#ifndef ZPP_SERIALIZER_FREESTANDING
        try {
            // Load the items.
            memory_view_input_archive::operator()(
                std::forward<Items>(items)...);
        } catch (...) {
            // Reset the offset back.
            reset(offset);
            throw;
        }
#else // ZPP_SERIALIZER_FREESTANDING

        // Load the items.
        if (auto result = memory_view_input_archive::operator()(
                std::forward<Items>(items)...);
            !result) {
            // Reset the offset back.
            reset(offset);
            return result;
        }
#endif

        // Erase the loaded elements.
        m_input->erase(m_input->begin(),
                       m_input->begin() + this->offset());

        // Reset to offset zero.
        reset();

#ifdef ZPP_SERIALIZER_FREESTANDING
        return freestanding::error{error::success};
#endif
    }

    /**
     * Returns the data pointer;
     */
    using base::data;

    /**
     * Returns the current offset in the data.
     */
    using base::offset;

    /**
     * Allow to reset offset for advanced use.
     */
    using base::reset;

private:
    /**
     * The input data.
     */
    std::vector<unsigned char> * m_input{};
};

#ifndef ZPP_SERIALIZER_FREESTANDING
/**
 * This class manages polymorphic type registration for serialization
 * process.
 */
template <typename Archive>
class registry
{
public:
    static_assert(!std::is_reference<Archive>::value,
                  "Disallows reference type for archive in registry");

    /**
     * Returns the global instance of the registry.
     */
    static registry & get_instance() noexcept
    {
        static registry registry;
        return registry;
    }

    /**
     * Add a serialization method for a given polymorphic type and id.
     */
    template <typename Type, id_type id>
    void add()
    {
        add<Type>(id);
    }

    /**
     * Adds a serialization method for a given polymorphic type and id.
     */
    template <typename Type>
    void add(id_type id)
    {
        add(id,
            typeid(Type).name(),
            make_serialization_method<Archive, Type>());
    }

    /**
     * Add a serialization method for a given polymorphic type information
     * string and id. The behavior is undefined if the type isn't derived
     * from polymorphic.
     */
    void add(id_type id,
             std::string type_information_string,
             serialization_method_t<Archive> serialization_method)
    {
        // Lock the serialization method maps for write access.
        std::lock_guard<shared_mutex> lock(m_shared_mutex);

        // Add the serialization id to serialization method mapping.
        m_serialization_id_to_method.emplace(
            id, std::move(serialization_method));

        // Add the type information to to serialization id mapping.
        m_type_information_to_serialization_id.emplace(
            std::move(type_information_string), id);
    }

    /**
     * Serialize a polymorphic type, in case of a loading (input) archive.
     */
    template <typename...,
              typename ArchiveType = Archive,
              typename = typename ArchiveType::loading>
    void serialize(Archive & archive,
                   std::unique_ptr<polymorphic> & object)
    {
        id_type id{};

        // Load the serialization id.
        archive(id);

        // Lock the serialization method maps for read access.
        std::shared_lock<shared_mutex> lock(m_shared_mutex);

        // Find the serialization method.
        auto serialization_id_to_method_pair =
            m_serialization_id_to_method.find(id);
        if (m_serialization_id_to_method.end() ==
            serialization_id_to_method_pair) {
            throw undeclared_polymorphic_type_error(
                "Undeclared polymorphic serialization type error.");
        }

        // Fetch the serialization method.
        auto serialization_method =
            serialization_id_to_method_pair->second;

        // Unlock the serialization method maps.
        lock.unlock();

        // Serialize (load) the given object.
        serialization_method(archive, object);
    }

    /**
     * Serialize a polymorphic type, in case of a saving (output) archive.
     */
    template <typename...,
              typename ArchiveType = Archive,
              typename = typename ArchiveType::saving>
    void serialize(Archive & archive, const polymorphic & object)
    {
        // Lock the serialization method maps for read access.
        std::shared_lock<shared_mutex> lock(m_shared_mutex);

        // Find the serialization id.
        auto type_information_to_serialization_id_pair =
            m_type_information_to_serialization_id.find(
                typeid(object).name());
        if (m_type_information_to_serialization_id.end() ==
            type_information_to_serialization_id_pair) {
            throw undeclared_polymorphic_type_error(
                "Undeclared polymorphic serialization type error.");
        }

        // Fetch the serialization id.
        auto id = type_information_to_serialization_id_pair->second;

        // Find the serialization method.
        auto serialization_id_to_method_pair =
            m_serialization_id_to_method.find(id);
        if (m_serialization_id_to_method.end() ==
            serialization_id_to_method_pair) {
            throw undeclared_polymorphic_type_error(
                "Undeclared polymorphic serialization type error.");
        }

        // Fetch the serialization method.
        auto serialization_method =
            serialization_id_to_method_pair->second;

        // Unlock the serialization method maps.
        lock.unlock();

        // Serialize (save) the serialization id.
        archive(id);

        // Serialize (save) the given object.
        serialization_method(archive, object);
    }

private:
    /**
     * Default constructor, defaulted.
     */
    registry() = default;

private:
    /**
     * The shared mutex that protects the maps below.
     */
    shared_mutex m_shared_mutex;

    /**
     * A map between serialization id to method.
     */
    std::unordered_map<id_type, serialization_method_t<Archive>>
        m_serialization_id_to_method;

    /**
     * A map between type information string to serialization id.
     */
    std::unordered_map<std::string, id_type>
        m_type_information_to_serialization_id;
};     // registry
#endif // ZPP_SERIALIZER_FREESTANDING

/**
 * Serialize resizable containers, operates on loading (input) archives.
 */
template <
    typename Archive,
    typename Container,
    typename SizeType = size_type,
    typename...,
    typename = decltype(std::declval<Container &>().size()),
    typename = decltype(std::declval<Container &>().begin()),
    typename = decltype(std::declval<Container &>().end()),
    typename = decltype(std::declval<Container &>().resize(std::size_t())),
    typename = std::enable_if_t<
        std::is_class<typename Container::value_type>::value ||
        !std::is_base_of<
            std::random_access_iterator_tag,
            typename std::iterator_traits<
                typename Container::iterator>::iterator_category>::value>,
    typename = typename Archive::loading,
    typename = void,
    typename = void,
    typename = void,
    typename = void>
auto serialize(Archive & archive, Container & container)
{
    SizeType size{};

    // Fetch the number of items to load.
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(size);
#else
    if (auto result = archive(size); !result) {
        return result;
    }
#endif

    // Resize the container to match the size.
    container.resize(size);

    // Serialize all the items.
    for (auto & item : container) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(item);
#else
        if (auto result = archive(item); !result) {
            return result;
        }
#endif
    }

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize containers, operates on saving (output) archives.
 */
template <
    typename Archive,
    typename Container,
    typename SizeType = size_type,
    typename...,
    typename = decltype(std::declval<Container &>().size()),
    typename = decltype(std::declval<Container &>().begin()),
    typename = decltype(std::declval<Container &>().end()),
    typename = std::enable_if_t<
        std::is_class<typename Container::value_type>::value ||
        !std::is_base_of<
            std::random_access_iterator_tag,
            typename std::iterator_traits<
                typename Container::iterator>::iterator_category>::value ||
        !detail::has_data_member_function<Container>::value>,
    typename = typename Archive::saving,
    typename = void,
    typename = void,
    typename = void,
    typename = void,
    typename = void>
auto serialize(Archive & archive, const Container & container)
{
#ifndef ZPP_SERIALIZER_FREESTANDING
    // Save the container size.
    archive(static_cast<SizeType>(container.size()));
#else
    if (auto result = archive(static_cast<SizeType>(container.size()));
        !result) {
        return result;
    }
#endif

    // Serialize all the items.
    for (auto & item : container) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(item);
#else
        if (auto result = archive(item); !result) {
            return result;
        }
#endif
    }

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize view containers, operates on loading (input) archives.
 */
template <
    typename Archive,
    typename Container,
    typename SizeType = size_type,
    typename...,
    typename = decltype(std::declval<Container &>().size()),
    typename = decltype(std::declval<Container &>().begin()),
    typename = decltype(std::declval<Container &>().end()),
    typename =
        std::enable_if_t<std::is_trivially_destructible<Container>::value>,
    typename = std::enable_if_t<
        std::is_class<typename Container::value_type>::value ||
        !std::is_base_of<
            std::random_access_iterator_tag,
            typename std::iterator_traits<
                typename Container::iterator>::iterator_category>::value>,
    typename = typename Archive::loading,
    typename = void,
    typename = void,
    typename = void,
    typename = void,
    typename = void>
auto serialize(Archive & archive, Container & container)
{
    SizeType size{};

    // Fetch the number of items to load.
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(size);
#else
    if (auto result = archive(size); !result) {
        return result;
    }
#endif

    // Check size.
    if (size > container.size()) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        throw out_of_range("View type container out of range.");
#else
        return freestanding::error{error::out_of_range};
#endif
    }

    // Resize the view container to match the size.
    container = {container.data(), size};

    // Serialize all the items.
    for (auto & item : container) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(item);
#else
        if (auto result = archive(item); !result) {
            return result;
        }
#endif

#ifdef ZPP_SERIALIZER_FREESTANDING
        return freestanding::error{error::success};
#endif
    }
}

/**
 * Serialize resizable, continuous containers, of fundamental or
 * enumeration types. Operates on loading (input) archives.
 */
template <
    typename Archive,
    typename Container,
    typename SizeType = size_type,
    typename...,
    typename = decltype(std::declval<Container &>().size()),
    typename = decltype(std::declval<Container &>().begin()),
    typename = decltype(std::declval<Container &>().end()),
    typename = decltype(std::declval<Container &>().resize(std::size_t())),
    typename = decltype(std::declval<Container &>().data()),
    typename = std::enable_if_t<
        std::is_fundamental<typename Container::value_type>::value ||
        std::is_enum<typename Container::value_type>::value>,
    typename = std::enable_if_t<std::is_base_of<
        std::random_access_iterator_tag,
        typename std::iterator_traits<
            typename Container::iterator>::iterator_category>::value>,
    typename = typename Archive::loading,
    typename = void,
    typename = void,
    typename = void,
    typename = void>
auto serialize(Archive & archive, Container & container)
{
    SizeType size{};

    // Fetch the number of items to load.
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(size);
#else
    if (auto result = archive(size); !result) {
        return result;
    }
#endif

    // Resize the container to match the size.
    container.resize(size);

    // If the size is zero, return.
    if (!size) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        return;
#else
        return freestanding::error{error::success};
#endif
    }

    // Serialize the bytes data.
    return archive(as_bytes(std::addressof(container[0]),
                            static_cast<SizeType>(container.size())));
}

/**
 * Serialize continuous containers, of fundamental or
 * enumeration types. Operates on saving (output) archives.
 */
template <typename Archive,
          typename Container,
          typename SizeType = size_type,
          typename...,
          typename = decltype(std::declval<Container &>().size()),
          typename = decltype(std::declval<Container &>().begin()),
          typename = decltype(std::declval<Container &>().end()),
          typename = decltype(std::declval<Container &>().data()),
          typename = std::enable_if_t<
              std::is_fundamental<typename Container::value_type>::value ||
              std::is_enum<typename Container::value_type>::value>,
          typename = std::enable_if_t<std::is_base_of<
              std::random_access_iterator_tag,
              typename std::iterator_traits<typename Container::iterator>::
                  iterator_category>::value>,
          typename = typename Archive::saving,
          typename = void,
          typename = void,
          typename = void,
          typename = void,
          typename = void>
auto serialize(Archive & archive, const Container & container)
{
    // The container size.
    auto size = static_cast<SizeType>(container.size());

    // Save the container size.
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(size);
#else
    if (auto result = archive(size); !result) {
        return result;
    }
#endif

    // If the size is zero, return.
    if (!size) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        return;
#else
        return freestanding::error{error::success};
#endif
    }

    // Serialize the bytes data.
    return archive(as_bytes(std::addressof(container[0]),
                            static_cast<SizeType>(container.size())));
}

/**
 * Serialize continuous view containers, of fundamental or
 * enumeration types. Operates on loading (input) archives.
 */
template <typename Archive,
          typename Container,
          typename SizeType = size_type,
          typename...,
          typename = decltype(std::declval<Container &>().size()),
          typename = decltype(std::declval<Container &>().begin()),
          typename = decltype(std::declval<Container &>().end()),
          typename = decltype(std::declval<Container &>().data()),
          typename = std::enable_if_t<
              std::is_trivially_destructible<Container>::value>,
          typename = std::enable_if_t<
              std::is_fundamental<typename Container::value_type>::value ||
              std::is_enum<typename Container::value_type>::value>,
          typename = std::enable_if_t<std::is_base_of<
              std::random_access_iterator_tag,
              typename std::iterator_traits<typename Container::iterator>::
                  iterator_category>::value>,
          typename = typename Archive::loading,
          typename = void,
          typename = void,
          typename = void,
          typename = void,
          typename = void>
auto serialize(Archive & archive, Container & container)
{
    SizeType size{};

    // Fetch the number of items to load.
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(size);
#else
    if (auto result = archive(size); !result) {
        return result;
    }
#endif

    // Check the size.
    if (size > container.size()) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        throw out_of_range("View type container out of range.");
#else
        return freestanding::error{error::out_of_range};
#endif
    }

    // Resize the view container to match the size.
    container = {container.data(), size};

    // If the size is zero, return.
    if (!size) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        return;
#else
        return freestanding::error{error::success};
#endif
    }

    // Serialize the bytes data.
    return archive(as_bytes(std::addressof(container[0]),
                            static_cast<SizeType>(container.size())));
}

/**
 * Serialize Associative and UnorderedAssociative containers, operates on
 * loading (input) archives.
 */
template <typename Archive,
          typename Container,
          typename SizeType = size_type,
          typename...,
          typename = decltype(std::declval<Container &>().size()),
          typename = decltype(std::declval<Container &>().begin()),
          typename = decltype(std::declval<Container &>().end()),
          typename = typename Container::value_type,
          typename = typename Container::key_type,
          typename = typename Archive::loading>
auto serialize(Archive & archive, Container & container)
{
    SizeType size{};

    // Fetch the number of items to load.
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(size);
#else
    if (auto result = archive(size); !result) {
        return result;
    }
#endif

    // Serialize all the items.
    for (SizeType i{}; i < size; ++i) {
        // Deduce the container item type.
        using item_type =
            detail::container_nonconst_value_type_t<Container>;

        // Create just enough storage properly aligned for one item.
        std::aligned_storage_t<sizeof(item_type), alignof(item_type)>
            storage;

        // Create the object at the storage.
        std::unique_ptr<item_type, void (*)(item_type *)> object(
            access::placement_new<item_type>(std::addressof(storage)),
            [](auto pointer) { access::destruct(*pointer); });

        // Serialize the object.
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(*object);
#else
        if (auto result = archive(*object); !result) {
            return result;
        }
#endif

        // Insert the item to the container.
        container.insert(std::move(*object));
    }

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize arrays, operates on loading (input) archives.
 * This overload is for non fundamental non enumeration types.
 */
template <typename Archive,
          typename Item,
          std::size_t size,
          typename...,
          typename = std::enable_if_t<!std::is_fundamental<Item>::value &&
                                      !std::is_enum<Item>::value>,
          typename = typename Archive::loading>
auto serialize(Archive & archive, Item (&array)[size])
{
    // Serialize every item.
    for (auto & item : array) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(item);
#else
        if (auto result = archive(item); !result) {
            return result;
        }
#endif
    }

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize arrays, operates on loading (input) archives.
 * This overload is for fundamental or enumeration types.
 */
template <typename Archive,
          typename Item,
          std::size_t size,
          typename...,
          typename = std::enable_if_t<std::is_fundamental<Item>::value ||
                                      std::is_enum<Item>::value>,
          typename = typename Archive::loading,
          typename = void>
auto serialize(Archive & archive, Item (&array)[size])
{
    return archive(as_bytes(array, size));
}

/**
 * Serialize arrays, operates on saving (output) archives.
 * This overload is for non fundamental non enumeration types.
 */
template <typename Archive,
          typename Item,
          std::size_t size,
          typename...,
          typename = std::enable_if_t<!std::is_fundamental<Item>::value &&
                                      !std::is_enum<Item>::value>,
          typename = typename Archive::saving>
auto serialize(Archive & archive, const Item (&array)[size])
{
    // Serialize every item.
    for (auto & item : array) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(item);
#else
        if (auto result = archive(item); !result) {
            return result;
        }
#endif
    }

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize arrays, operates on saving (output) archives.
 * This overload is for fundamental or enumeration types.
 */
template <typename Archive,
          typename Item,
          std::size_t size,
          typename...,
          typename = std::enable_if_t<std::is_fundamental<Item>::value ||
                                      std::is_enum<Item>::value>,
          typename = typename Archive::saving,
    typename = void>
auto serialize(Archive & archive, const Item (&array)[size])
{
    return archive(as_bytes(array, size));
}

/**
 * Serialize std::array, operates on loading (input) archives.
 * This overload is for non fundamental non enumeration types.
 */
template <typename Archive,
          typename Item,
          std::size_t size,
          typename...,
          typename = std::enable_if_t<!std::is_fundamental<Item>::value &&
                                      !std::is_enum<Item>::value>,
          typename = typename Archive::loading>
auto serialize(Archive & archive, std::array<Item, size> & array)
{
    // Serialize every item.
    for (auto & item : array) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(item);
#else
        if (auto result = archive(item); !result) {
            return result;
        }
#endif
    }

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize std::array, operates on loading (input) archives.
 * This overload is for fundamental or enumeration types.
 */
template <typename Archive,
          typename Item,
          std::size_t size,
          typename...,
          typename = std::enable_if_t<std::is_fundamental<Item>::value ||
                                      std::is_enum<Item>::value>,
          typename = typename Archive::loading,
          typename = void>
auto serialize(Archive & archive, std::array<Item, size> & array)
{
    return archive(as_bytes(std::addressof(array[0]), size));
}

/**
 * Serialize std::array, operates on saving (output) archives.
 * This overload is for non fundamental non enumeration types.
 */
template <typename Archive,
          typename Item,
          std::size_t size,
          typename...,
          typename = std::enable_if_t<!std::is_fundamental<Item>::value &&
                                      !std::is_enum<Item>::value>,
          typename = typename Archive::saving>
auto serialize(Archive & archive, const std::array<Item, size> & array)
{
    // Serialize every item.
    for (auto & item : array) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(item);
#else
        if (auto result = archive(item); !result) {
            return result;
        }
#endif
    }

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize std::array, operates on saving (output) archives.
 * This overload is for fundamental or enumeration types.
 */
template <typename Archive,
          typename Item,
          std::size_t size,
          typename...,
          typename = std::enable_if_t<std::is_fundamental<Item>::value ||
                                      std::is_enum<Item>::value>,
          typename = typename Archive::saving,
          typename = void>
auto serialize(Archive & archive, const std::array<Item, size> & array)
{
    return archive(as_bytes(std::addressof(array[0]), size));
}

/**
 * Serialize std::pair, operates on loading (input) archives.
 */
template <typename Archive,
          typename First,
          typename Second,
          typename...,
          typename = typename Archive::loading>
auto serialize(Archive & archive, std::pair<First, Second> & pair)
{
    // Serialize first, then second.
    return archive(pair.first, pair.second);
}

/**
 * Serialize std::pair, operates on saving (output) archives.
 */
template <typename Archive,
          typename First,
          typename Second,
          typename...,
          typename = typename Archive::saving>
auto serialize(Archive & archive, const std::pair<First, Second> & pair)
{
    // Serialize first, then second.
    return archive(pair.first, pair.second);
}

/**
 * Serialize std::tuple, operates on loading (input) archives.
 */
template <typename Archive,
          typename... TupleItems,
          typename = typename Archive::loading>
auto serialize(Archive & archive, std::tuple<TupleItems...> & tuple)
{
    // Delegate to a helper function with an index sequence.
    return serialize(
        archive, tuple, std::make_index_sequence<sizeof...(TupleItems)>());
}

/**
 * Serialize std::tuple, operates on saving (output) archives.
 */
template <typename Archive,
          typename... TupleItems,
          typename = typename Archive::saving>
auto serialize(Archive & archive, const std::tuple<TupleItems...> & tuple)
{
    // Delegate to a helper function with an index sequence.
    return serialize(
        archive, tuple, std::make_index_sequence<sizeof...(TupleItems)>());
}

/**
 * Serialize std::tuple, operates on loading (input) archives.
 * This overload serves as a helper function that accepts an index
 * sequence.
 */
template <typename Archive,
          typename... TupleItems,
          std::size_t... Indices,
          typename = typename Archive::loading>
auto serialize(Archive & archive,
               std::tuple<TupleItems...> & tuple,
               std::index_sequence<Indices...>)
{
    return archive(std::get<Indices>(tuple)...);
}

/**
 * Serialize std::tuple, operates on saving (output) archives.
 * This overload serves as a helper function that accepts an index
 * sequence.
 */
template <typename Archive,
          typename... TupleItems,
          std::size_t... Indices,
          typename = typename Archive::saving>
auto serialize(Archive & archive,
               const std::tuple<TupleItems...> & tuple,
               std::index_sequence<Indices...>)
{
    return archive(std::get<Indices>(tuple)...);
}

#if __cplusplus >= 201703L
/**
 * Serialize std::optional, operates on loading (input) archives.
 */
template <typename Archive,
          typename Type,
          typename...,
          typename = typename Archive::loading>
auto serialize(Archive & archive, std::optional<Type> & optional)
{
    // Load whether has value.
    bool has_value{};
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(has_value);
#else
    if (auto result = archive(has_value); !result) {
        return result;
    }
#endif

    // If does not have a value.
    if (!has_value) {
        optional = std::nullopt;
#ifndef ZPP_SERIALIZER_FREESTANDING
        return;
#else
        return freestanding::error{error::success};
#endif
    }

    // If the type is default constructible.
    if constexpr (std::is_default_constructible_v<Type>) {
        // Create the value if does not exist.
        if (!optional) {
            optional = Type{};
        }

        // Load the value.
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(*optional);
#else
        if (auto result = archive(*optional); !result) {
            return result;
        }
#endif
    } else {
        // The object storage.
        std::aligned_storage_t<sizeof(Type), alignof(Type)> storage;

        // Create the object at the storage.
        std::unique_ptr<Type, void (*)(Type *)> object(
            access::placement_new<Type>(std::addressof(storage)),
            [](auto pointer) { access::destruct(*pointer); });

        // Load the object.
#ifndef ZPP_SERIALIZER_FREESTANDING
        archive(*object);
#else
        if (auto result = archive(*optional); !result) {
            return result;
        }
#endif

        // Assign the loaded object.
        optional = std::move(*object);
    }

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize std::optional, operates on saving (output) archives.
 */
template <typename Archive,
          typename Type,
          typename...,
          typename = typename Archive::saving>
auto serialize(Archive & archive, const std::optional<Type> & optional)
{
    // Save has value.
    bool has_value = optional.has_value();

    // If has value, save it.
    if (has_value) {
        return archive(has_value, *optional);
    } else {
        return archive(has_value);
    }
}

/**
 * Serialize std::variant, operates on loading (input) archives.
 */
template <typename Archive,
          typename... Types,
          typename = typename Archive::loading>
auto serialize(Archive & archive, std::variant<Types...> & variant)
{
    // Test for maximum number of types.
    static_assert(sizeof...(Types) < 0xff, "Max variant types reached.");

    // The variant index.
    unsigned char index{};

    // Load the index.
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(index);
#else
    if (auto result = archive(index); !result) {
        return result;
    }
#endif

    // Check that loaded index is inside bounds.
    if (index >= sizeof...(Types)) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        throw variant_index_out_of_range("Variant index out of range");
#else
        return freestanding::error{error::out_of_range};
#endif
    }

    // The variant type.
    using variant_type = std::variant<Types...>;

    // Loader type.
    using loader_type =
        void (*)(Archive & archive, variant_type & variant);

    // Loaders per variant index.
    static constexpr loader_type loaders[] = {[](auto & archive,
                                                 auto & variant) {
        // If the type is default constructible.
        if constexpr (std::is_default_constructible_v<Types>) {
            // If does not have the needed type, assign it.
            if (!std::get_if<Types>(&variant)) {
                variant = Types{};
            }

            // Load the value.
            return archive(*std::get_if<Types>(&variant));
        } else {
            // The object storage.
            std::aligned_storage_t<sizeof(Types), alignof(Types)> storage;

            // Create the object at the storage.
            std::unique_ptr<Types, void (*)(Types *)> object(
                access::placement_new<Types>(std::addressof(storage)),
                [](auto pointer) { access::destruct(*pointer); });

            // Load the object.
#ifndef ZPP_SERIALIZER_FREESTANDING
            archive(*object);
#else
            if (auto result = archive(*object); !result) {
                return result;
            }
#endif

            // Assign the loaded object.
            variant = std::move(*object);
        }
    }...};

    // Execute the appropriate loader.
    return loaders[index](archive, variant);
}

/**
 * Serialize std::variant, operates on saving (output) archives.
 */
template <typename Archive,
          typename... Types,
          typename = typename Archive::saving>
auto serialize(Archive & archive, const std::variant<Types...> & variant)
{
    // Test for maximum number of types.
    static_assert(sizeof...(Types) < 0xff, "Max variant types reached.");

    // The variant index.
    auto variant_index = variant.index();

    // Disallow serializations of valueless variant.
    if (std::variant_npos == variant_index) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        throw attempt_to_serialize_valueless_variant(
            "Cannot serialize a valueless variant.");
#else
        return freestanding::error{error::variant_is_valueless};
#endif
    }

    // The index to save.
    auto index = static_cast<unsigned char>(variant_index & 0xff);

    // Save the variant object.
    return std::visit(
        [index, &archive](auto & object) {
            return archive(index, object);
        },
        variant);
}
#endif

/**
 * Serialize std::unique_ptr of non polymorphic, in case of a loading
 * (input) archive.
 */
template <typename Archive,
          typename Type,
          typename...,
          typename =
              std::enable_if_t<!std::is_base_of<polymorphic, Type>::value>,
          typename = typename Archive::loading>
auto serialize(Archive & archive, std::unique_ptr<Type> & object)
{
    // Construct a new object.
    auto loaded_object = access::make_unique<Type>();

    // Serialize the object.
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(*loaded_object);
#else
    if (auto result = archive(*loaded_object); !result) {
        return result;
    }
#endif

    // Transfer the object.
    object.reset(loaded_object.release());

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize std::unique_ptr of non polymorphic, in case of a saving
 * (output) archive.
 */
template <typename Archive,
          typename Type,
          typename...,
          typename =
              std::enable_if_t<!std::is_base_of<polymorphic, Type>::value>,
          typename = typename Archive::saving>
auto serialize(Archive & archive, const std::unique_ptr<Type> & object)
{
    // Prevent serialization of null pointers.
    if (nullptr == object) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        throw attempt_to_serialize_null_pointer_error(
            "Attempt to serialize null pointer.");
#else
        return freestanding::error{error::null_pointer_serialization};
#endif
    }

    // Serialize the object.
    return archive(*object);
}

#ifndef ZPP_SERIALIZER_FREESTANDING
/**
 * Serialize std::unique_ptr of polymorphic, in case of a loading (input)
 * archive.
 */
template <
    typename Archive,
    typename Type,
    typename...,
    typename = std::enable_if_t<std::is_base_of<polymorphic, Type>::value>,
    typename = typename Archive::loading,
    typename = void>
void serialize(Archive & archive, std::unique_ptr<Type> & object)
{
    std::unique_ptr<polymorphic> loaded_type;

    // Get the instance of the polymorphic registry.
    auto & registry_instance = registry<Archive>::get_instance();

    // Serialize the object using the registry.
    registry_instance.serialize(archive, loaded_type);

    try {
        // Check if the loaded type is convertible to Type.
        object.reset(&dynamic_cast<Type &>(*loaded_type));

        // Release the object.
        loaded_type.release();
    } catch (const std::bad_cast &) {
        // The loaded type was not convertible to Type.
        throw polymorphic_type_mismatch_error(
            "Polymorphic serialization type mismatch.");
    }
}

/**
 * Serialize std::unique_ptr of polymorphic, in case of a saving (output)
 * archive.
 */
template <
    typename Archive,
    typename Type,
    typename...,
    typename = std::enable_if_t<std::is_base_of<polymorphic, Type>::value>,
    typename = typename Archive::saving,
    typename = void>
void serialize(Archive & archive, const std::unique_ptr<Type> & object)
{
    // Prevent serialization of null pointers.
    if (nullptr == object) {
        throw attempt_to_serialize_null_pointer_error(
            "Attempt to serialize null pointer.");
    }

    // Get the instance of the polymorphic registry.
    auto & registry_instance = registry<Archive>::get_instance();

    // Serialize the object using the registry.
    registry_instance.serialize(archive, *object);
}
#endif // ZPP_SERIALIZER_FREESTANDING

/**
 * Serialize std::shared_ptr of non polymorphic, in case of a loading
 * (input) archive.
 */
template <typename Archive,
          typename Type,
          typename...,
          typename =
              std::enable_if_t<!std::is_base_of<polymorphic, Type>::value>,
          typename = typename Archive::loading>
auto serialize(Archive & archive, std::shared_ptr<Type> & object)
{
    // Construct a new object.
    auto loaded_object = access::make_unique<Type>();

    // Serialize the object.
#ifndef ZPP_SERIALIZER_FREESTANDING
    archive(*loaded_object);
#else
    if (auto result = archive(*loaded_object); !result) {
        return result;
    }
#endif

    // Transfer the object.
    object.reset(loaded_object.release());

#ifdef ZPP_SERIALIZER_FREESTANDING
    return freestanding::error{error::success};
#endif
}

/**
 * Serialize std::shared_ptr of non polymorphic, in case of a saving
 * (output) archive.
 */
template <typename Archive,
          typename Type,
          typename...,
          typename =
              std::enable_if_t<!std::is_base_of<polymorphic, Type>::value>,
          typename = typename Archive::saving>
auto serialize(Archive & archive, const std::shared_ptr<Type> & object)
{
    // Prevent serialization of null pointers.
    if (nullptr == object) {
#ifndef ZPP_SERIALIZER_FREESTANDING
        throw attempt_to_serialize_null_pointer_error(
            "Attempt to serialize null pointer.");
#else
        return freestanding::error{error::null_pointer_serialization};
#endif
    }

    // Serialize the object.
    return archive(*object);
}

/**
 * Represents a container object with specific
 * size type requirements.
 */
template <typename SizeType, typename Container>
class sized_container
{
public:
    /**
     * Must be class type.
     */
    static_assert(std::is_class<Container>::value,
                  "Container must be a class type.");

    /**
     * Must be unsigned integral type.
     */
    static_assert(std::is_unsigned<SizeType>::value,
                  "Size must be an unsigned integral type.");

    /*
     * Construct the sized container.
     */
    explicit sized_container(Container & container) : container(container)
    {
    }

    /**
     * Call serialize directly with the size type parameter.
     */
    template <typename Archive, typename Self>
    static auto serialize(Archive & archive, Self & self)
    {
        using zpp::serializer::serialize;
        return serialize<Archive, Container, SizeType>(archive,
                                                       self.container);
    }

    /**
     * The wrapped container type.
     */
    Container & container;
};

/**
 * Creates a wrapper object of sized_container to
 * allow serialization with specific size type requirements.
 */
template <typename SizeType, typename Container>
auto size_is(Container && container)
{
    return sized_container<SizeType, std::remove_reference_t<Container>>(
        container);
}

#ifndef ZPP_SERIALIZER_FREESTANDING
/**
 * Serialize std::shared_ptr of polymorphic, in case of a loading (input)
 * archive.
 */
template <
    typename Archive,
    typename Type,
    typename...,
    typename = std::enable_if_t<std::is_base_of<polymorphic, Type>::value>,
    typename = typename Archive::loading,
    typename = void>
void serialize(Archive & archive, std::shared_ptr<Type> & object)
{
    std::unique_ptr<polymorphic> loaded_type;

    // Get the instance of the polymorphic registry.
    auto & registry_instance = registry<Archive>::get_instance();

    // Serialize the object using the registry.
    registry_instance.serialize(archive, loaded_type);

    try {
        // Check if the loaded type is convertible to Type.
        object.reset(&dynamic_cast<Type &>(*loaded_type));

        // Release the object.
        loaded_type.release();
    } catch (const std::bad_cast &) {
        // The loaded type was not convertible to Type.
        throw polymorphic_type_mismatch_error(
            "Polymorphic serialization type mismatch.");
    }
}

/**
 * Serialize std::shared_ptr of polymorphic, in case of a saving (output)
 * archive.
 */
template <
    typename Archive,
    typename Type,
    typename...,
    typename = std::enable_if_t<std::is_base_of<polymorphic, Type>::value>,
    typename = typename Archive::saving,
    typename = void>
void serialize(Archive & archive, const std::shared_ptr<Type> & object)
{
    // Prevent serialization of null pointers.
    if (nullptr == object) {
        throw attempt_to_serialize_null_pointer_error(
            "Attempt to serialize null pointer.");
    }

    // Get the instance of the polymorphic registry.
    auto & registry_instance = registry<Archive>::get_instance();

    // Serialize the object using the registry.
    registry_instance.serialize(archive, *object);
}

/**
 * Serialize types wrapped with polymorphic_wrapper,
 * which is supported only for saving (output) archives.
 * Usually used with the as_polymorphic facility.
 */
template <typename Archive,
          typename Type,
          typename...,
          typename = typename Archive::saving>
void serialize(Archive & archive, const polymorphic_wrapper<Type> & object)
{
    // Get the instance of the polymorphic registry.
    auto & registry_instance = registry<Archive>::get_instance();

    // Serialize using the registry.
    registry_instance.serialize(archive, *object);
}

/**
 * A meta container that holds a sequence of archives.
 */
template <typename... Archives>
struct archive_sequence
{
};

/**
 * The built in archives.
 */
using builtin_archives = archive_sequence<memory_view_input_archive,
                                          basic_memory_output_archive>;

/**
 * Makes a meta pair of type and id.
 */
template <typename Type, id_type id>
struct make_type;

/**
 * Registers user defined polymorphic types to serialization registry.
 */
template <typename... ExtraTypes>
class register_types;

/**
 * A no operation class, registers an empty list of types.
 */
template <>
class register_types<>
{
};

/**
 * Registers user defined polymorphic types to serialization registry.
 */
template <typename Type, id_type id, typename... ExtraTypes>
class register_types<make_type<Type, id>, ExtraTypes...>
    : private register_types<ExtraTypes...>
{
public:
    /**
     * Registers the type to the built in archives of the serializer.
     */
    register_types() noexcept : register_types(builtin_archives())
    {
    }

    /**
     * Registers the type to every archive in the given archive sequence.
     */
    template <typename... Archives>
    register_types(archive_sequence<Archives...> archives) noexcept
    {
        register_type_to_archives(archives);
    }

private:
    /**
     * Registers the type to every archive in the given archive sequence.
     */
    template <typename Archive, typename... Archives>
    void register_type_to_archives(
        archive_sequence<Archive, Archives...>) noexcept
    {
        // Register the type to the first archive.
        register_type_to_archive<Archive>();

        // Register the type to the other archives.
        register_type_to_archives(archive_sequence<Archives...>());
    }

    /**
     * Registers the type to an empty archive sequence - does nothing.
     */
    void register_type_to_archives(archive_sequence<>) noexcept
    {
    }

    /**
     * Registers the type to the given archive.
     * Must throw no exceptions since this will most likely execute
     * during static construction.
     * The effect of failure is that the type will not be registered,
     * it will be detected during runtime.
     */
    template <typename Archive>
    void register_type_to_archive() noexcept
    {
        try {
            registry<Archive>::get_instance().template add<Type, id>();
        } catch (...) {
        }
    }
}; // register_types

/**
 * Accepts a name and returns its serialization id.
 * We return the first 8 bytes of the sha1 on the given name.
 */
template <std::size_t size>
constexpr id_type make_id(const char (&name)[size])
{
    // Initialize constants.
    std::uint32_t h0 = 0x67452301u;
    std::uint32_t h1 = 0xEFCDAB89u;
    std::uint32_t h2 = 0x98BADCFEu;
    std::uint32_t h3 = 0x10325476u;
    std::uint32_t h4 = 0xC3D2E1F0u;

    // Initialize the message size in bits.
    std::uint64_t message_size = (size - 1) * 8;

    // Calculate the size aligned to 64 bytes (512 bits).
    constexpr std::size_t aligned_message_size =
        (((size + sizeof(std::uint64_t)) + 63) / 64) * 64;

    // Construct the pre-processed message.
    std::uint32_t preprocessed_message[aligned_message_size /
                                       sizeof(std::uint32_t)] = {};
    for (std::size_t i{}; i < size - 1; ++i) {
        preprocessed_message[i / 4] |= detail::swap_byte_order(
            std::uint32_t(name[i])
            << ((sizeof(std::uint32_t) - 1 - (i % 4)) * 8));
    }

    // Append the byte 0x80.
    preprocessed_message[(size - 1) / 4] |= detail::swap_byte_order(
        std::uint32_t(0x80)
        << ((sizeof(std::uint32_t) - 1 - ((size - 1) % 4)) * 8));

    // Append the length in bits, in 64 bit, big endian.
    preprocessed_message[(aligned_message_size / sizeof(std::uint32_t)) -
                         2] =
        detail::swap_byte_order(std::uint32_t(message_size >> 32));
    preprocessed_message[(aligned_message_size / sizeof(std::uint32_t)) -
                         1] =
        detail::swap_byte_order(std::uint32_t(message_size));

    // Process the message in successive 512-bit chunks.
    for (std::size_t i{};
         i < (aligned_message_size / sizeof(std::uint32_t));
         i += 16) {
        std::uint32_t w[80] = {};

        // Set the value of w.
        for (std::size_t j{}; j < 16; ++j) {
            w[j] = preprocessed_message[i + j];
        }

        // Extend the sixteen 32-bit words into eighty 32-bit words.
        for (std::size_t j = 16; j < 80; ++j) {
            w[j] = detail::swap_byte_order(detail::rotate_left(
                detail::swap_byte_order(w[j - 3] ^ w[j - 8] ^ w[j - 14] ^
                                        w[j - 16]),
                1));
        }

        // Initialize hash values for this chunk.
        auto a = h0;
        auto b = h1;
        auto c = h2;
        auto d = h3;
        auto e = h4;

        // Main loop.
        for (std::size_t j{}; j < 80; ++j) {
            std::uint32_t f{};
            std::uint32_t k{};
            if (j <= 19) {
                f = (b & c) | ((~b) & d);
                k = 0x5A827999u;
            } else if (j <= 39) {
                f = b ^ c ^ d;
                k = 0x6ED9EBA1u;
            } else if (j <= 59) {
                f = (b & c) | (b & d) | (c & d);
                k = 0x8F1BBCDCu;
            } else {
                f = b ^ c ^ d;
                k = 0xCA62C1D6u;
            }

            auto temp = detail::rotate_left(a, 5) + f + e + k +
                        detail::swap_byte_order(w[j]);
            e = d;
            d = c;
            c = detail::rotate_left(b, 30);
            b = a;
            a = temp;
        }

        // Add this chunk's hash to result so far.
        h0 += a;
        h1 += b;
        h2 += c;
        h3 += d;
        h4 += e;
    }

    // Produce the first 8 bytes of the hash in little endian.
    return detail::swap_byte_order((std::uint64_t(h0) << 32) | h1);
} // make_id
#endif // ZPP_SERIALIZER_FREESTANDING

} // namespace serializer
} // namespace zpp

#endif // ZPP_SERIALIZER_H
