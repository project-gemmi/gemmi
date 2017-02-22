
#include "cif.hh"
#include <cassert>
#include <iostream>

// based on tao/json/internal/escape.hh
void escape(std::ostream& os, const std::string& s) {
  static const char* h = "0123456789abcdef";
  const char* p = s.data();
  const char* l = p;
  const char* const e = s.data() + s.size();
  while (p != e) {
    const unsigned char c = *p;
    if (c == '\\') {
      os.write(l, p - l);
      l = ++p;
      os << "\\\\";
    } else if (c == '"') {
      os.write(l, p - l);
      l = ++p;
      os << "\\\"";
    } else if (c < 32) {
      os.write(l, p - l);
      l = ++p;
      switch ( c ) {
        case '\b': os << "\\b"; break;
        case '\f': os << "\\f"; break;
        case '\n': os << "\\n"; break;
        case '\r': os << "\\r"; break;
        case '\t': os << "\\t"; break;
        default: os << "\\u00" << h[(c & 0xf0) >> 4] << h[c & 0x0f];
      }
    } else if (c == 127) {
      os.write(l, p - l);
      l = ++p;
      os << "\\u007f";
    } else {
      ++p;
    }
  }
  os.write(l, p - l);
}

class JsonWriter {
public:
  JsonWriter(std::ostream& os) : os_(os), linesep_("\n") {}
  void write_json(const cif::Document& d);

private:
  std::ostream& os_;
  std::string linesep_;

  void write_string(const std::string& s) {
    os_.put('"');
    escape(os_, s);
    os_.put('"');
  }

  void write_item(const cif::Item& item) {
    switch (item.type) {
      case cif::ItemType::Value:
        write_string(item.tv.tag);
        os_ << ": ";
        write_string(cif::as_string(item.tv.value));
        break;
      case cif::ItemType::Loop: {
        size_t ncol = item.loop.tags.size();
        const auto& vals = item.loop.values;
        for (size_t i = 0; i < ncol; i++) {
          if (i != 0)
            os_ << "," << linesep_;
          write_string(item.loop.tags[i].tag);
          os_ << ": [";
          for (size_t j = i; j < vals.size(); j += ncol) {
            if (j != i)
              os_.put(',');
            write_string(cif::as_string(vals[j]));
          }
          os_.put(']');
        }
        break;
      }
      case cif::ItemType::Frame:
        write_map(item.frame.tag, item.frame.items);
        break;
    }
  }

  void write_map(const std::string& name, const std::vector<cif::Item>& items) {
    os_ << linesep_;
    write_string(name);
    size_t n = linesep_.size();
    os_ << ": {";
    linesep_.resize(n + 1, ' ');
    for (const cif::Item& item : items) {
      if (&item != &items[0])
        os_.put(',');
      os_ << linesep_;
      write_item(item);
    }
    os_ << linesep_ << "}";
    linesep_.resize(n);
  }
};

void JsonWriter::write_json(const cif::Document& d) {
  os_.put('{');
  for (const cif::Block& block : d.blocks) {
    if (&block != &d.blocks[0])
      os_.put(',');
    write_map(block.name, block.items);
  }
  os_ << "\n}\n";
}

int main(int argc, char **argv) {
  if (argc < 2)
    return 1;
  const char* filename = argv[1];
  cif::Document d;
  try {
    d.parse_file(filename);
  } catch (pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  JsonWriter(std::cout).write_json(d);
  return 0;
}

// vim:sw=2:ts=2:et
