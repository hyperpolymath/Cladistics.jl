# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# config.ncl — Nickel configuration for byte detection and validation
#
# This file provides structured configuration for byte-level corruption detection,
# matching the patterns defined in stdlib/ByteDetector.affine.
#
# Language: Nickel (NCL) - https://nickel-lang.org/

{
  # Byte detection configuration
  byte_detection = {
    # Version of this configuration schema
    version = "1.0.0",

    # BLOCKING patterns - cause gate failure
    blocking = {
      # C0 control characters (ISO C0: 0x00-0x1F)
      # Excludes: 0x09 (tab), 0x0A (LF), 0x0D (CR)
      c0_control = {
        description = "C0 control characters and NUL bytes - file corruption",

        # PCRE pattern for detection
        pattern = "\\x00|[\\x01-\\x08\\x0B\\x0C\\x0E-\\x1F]",

        # Detailed character ranges
        ranges = {
          nul = {
            hex = "0x00",
            char = "NUL",
            description = "Null byte - file corruption indicator"
          },
          soh_to_bs = {
            hex_range = "0x01-0x08",
            chars = ["SOH", "STX", "ETX", "EOT", "ENQ", "ACK", "BEL", "BS"],
            description = "Start of heading through backspace"
          },
          vt = {
            hex = "0x0B",
            char = "VT",
            description = "Vertical tab"
          },
          ff = {
            hex = "0x0C",
            char = "FF",
            description = "Form feed"
          },
          so_to_us = {
            hex_range = "0x0E-0x1F",
            chars = ["SO", "SI", "DLE", "DC1", "DC2", "DC3", "DC4", "NAK", "SYN", "ETB", "CAN", "EM", "SUB", "ESC", "FS", "GS", "RS", "US"],
            description = "Shift out through unit separator"
          }
        },

        # Enforcement
        enforcement = {
          action = "block",
          exit_code = 1,
          severity = "error",
          message = "C0 control characters or NUL bytes - file corruption, blocks the gate"
        },

        # Rationale
        rationale = "These characters indicate file corruption or binary data in text files. They are NEVER legitimate in source code."
      }
    },

    # ADVISORY patterns - warning only
    advisory = {
      # Invisible Unicode characters
      invisible_unicode = {
        description = "Invisible Unicode characters (typography, legitimate in prose)",

        # PCRE pattern with UTF-8 mode
        # Combines C0 controls (for completeness) with invisible Unicode
        pattern = "(*UTF)[\\x00-\\x08\\x0B\\x0C\\x0E-\\x1F\\x{a0}\\x{ad}\\x{200b}-\\x{200f}\\x{202a}-\\x{202f}\\x{2060}\\x{2066}-\\x{2069}\\x{feff}]",

        # Character categories
        categories = {
          nbsp = {
            unicode = "U+00A0",
            hex = "\\x{a0}",
            name = "Non-Breaking Space",
            legitimate_use = "Typography in documentation"
          },
          soft_hyphen = {
            unicode = "U+00AD",
            hex = "\\x{ad}",
            name = "Soft Hyphen",
            legitimate_use = "Word breaking in documentation"
          },
          zero_width_spaces = {
            unicode_range = "U+200B-U+200F",
            hex = "\\x{200b}-\\x{200f}",
            names = ["ZWSP", "ZWNJ", "ZWJ", "LRM", "RLM"],
            legitimate_use = "Complex text layout, RTL languages"
          },
          bidi_formatting = {
            unicode_range = "U+202A-U+202F",
            hex = "\\x{202a}-\\x{202f}",
            names = ["LRE", "RLE", "PDF", "LRO", "RLO", "NNBSP"],
            legitimate_use = "Bidirectional text formatting"
          },
          word_joiner = {
            unicode = "U+2060",
            hex = "\\x{2060}",
            name = "Word Joiner",
            legitimate_use = "Typography"
          },
          bidi_isolates = {
            unicode_range = "U+2066-U+2069",
            hex = "\\x{2066}-\\x{2069}",
            names = ["LRI", "RLI", "FSI", "PDI"],
            legitimate_use = "Modern bidirectional text"
          },
          bom = {
            unicode = "U+FEFF",
            hex = "\\x{feff}",
            name = "Byte Order Mark (BOM) / Zero Width No-Break Space",
            legitimate_use = "Legacy: zero-width no-break space (deprecated)",
            note = "In modern Unicode, ZWNBSP is deprecated; U+2060 (WJ) should be used instead"
          }
        },

        # Enforcement
        enforcement = {
          action = "warn",
          exit_code = 0,
          severity = "warning",
          message = "Invisible Unicode characters detected (zero-width space, BOM, NBSP, etc.)"
        },

        # Rationale
        rationale = "About 2,100 estate files carry invisible Unicode as legitimate typography in prose (AsciiDoc/Markdown). Advisory only, not blocking."
      },

      # Leading BOM check (separate from general BOM detection)
      leading_bom = {
        description = "Byte Order Mark at file start (should use UTF-8 without BOM)",

        # Hex patterns for different BOM types
        patterns = {
          utf8 = {
            name = "UTF-8 BOM",
            hex_bytes = "EF BB BF",
            pattern = "ef bb bf"
          },
          utf16_le = {
            name = "UTF-16 Little Endian BOM",
            hex_bytes = "FF FE",
            pattern = "ff fe"
          },
          utf16_be = {
            name = "UTF-16 Big Endian BOM",
            hex_bytes = "FE FF",
            pattern = "fe ff"
          },
          utf32_le = {
            name = "UTF-32 Little Endian BOM",
            hex_bytes = "FF FE 00 00",
            pattern = "ff fe 00 00"
          },
          utf32_be = {
            name = "UTF-32 Big Endian BOM",
            hex_bytes = "00 00 FE FF",
            pattern = "00 00 fe ff"
          }
        },

        # Combined regex for detection
        combined_pattern = "ef bb bf|ff fe|fe ff|00 00 fe ff|ff fe 00 00",

        # Detection method
        detection = {
          command = "head -c 3 FILE | od -An -tx1 | grep -qE '^ *(ef bb bf|ff fe|fe ff|00 00 fe ff|ff fe 00 00)'",
          description = "Check first 3-4 bytes of file for BOM"
        },

        # Enforcement
        enforcement = {
          action = "warn",
          exit_code = 0,
          severity = "warning",
          message = "Leading BOM (Byte Order Mark) detected - should use UTF-8 without BOM"
        },

        # Rationale
        rationale = "Source files should use UTF-8 without BOM. BOM causes issues with parsers, compilers, and version control. Unix convention is no BOM."
      }
    },

    # File types to scan
    file_types = [
      "*.rs",     # Rust
      "*.ex",     # Elixir
      "*.exs",    # Elixir scripts
      "*.res",    # ReScript
      "*.js",     # JavaScript
      "*.ts",     # TypeScript
      "*.json",   # JSON
      "*.toml",   # TOML
      "*.yml",    # YAML
      "*.yaml",   # YAML
      "*.md",     # Markdown
      "*.adoc",   # AsciiDoc
      "*.idr",    # Idris
      "*.zig",    # Zig
      "*.v",      # V / Coq
      "*.jl",     # Julia
      "*.gleam",  # Gleam
      "*.hs",     # Haskell
      "*.ml",     # OCaml
      "*.sh",     # Shell scripts
      "*.affine", # Affine patterns
      "*.ncl"     # Nickel config
    ],

    # Directories to exclude from scanning
    exclude_paths = [
      ".git",
      "node_modules",
      ".deno",
      "target",
      "_build",
      "deps",
      "external_corpora",
      ".lake"
    ],

    # Grep options
    grep_options = {
      pcre_mode = true,           # Use -P for PCRE
      text_mode = true,           # Use -a to treat binary as text
      recursive = true,           # Use -r for recursive
      files_with_matches = true,  # Use -l for file names only
      quiet = true                # Use -q for exit code only (in blocking check)
    }
  },

  # Enforcement policy
  enforcement_policy = {
    version = "1.0.0",
    owner_ruling = "2026-08-28",
    owner = "Jonathan D.A. Jewell",

    summary = {
      blocking = "C0 control characters and NUL bytes only",
      advisory = "Invisible Unicode (NBSP, BOM, zero-width spaces)",
      rationale = "About 2,100 estate files carry legitimate invisible Unicode in prose. Only C0/NUL indicates corruption."
    }
  },

  # Tool compatibility
  compatibility = {
    pcre = {
      version = "PCRE 8.0+",
      required_features = ["UTF-8 mode", "hex escapes"],
      notes = "(*UTF) prefix for UTF-8 mode is PCRE-specific"
    },
    grep = {
      version = "GNU grep 2.5+",
      required_flags = ["-P", "-a"],
      notes = "Requires PCRE support (grep -P)"
    },
    locale = {
      independence = true,
      notes = "Uses hex escapes instead of character classes for locale independence"
    }
  },

  # References
  references = {
    c0_standard = "ISO/IEC 6429:1992 (C0 control characters)",
    unicode_standard = "Unicode Standard, Chapter 16 (Special Areas and Format Characters)",
    utf8_rfc = "RFC 3629 (UTF-8, Section 6: Byte Order Mark)",
    posix = "POSIX.1-2017 (grep, regular expressions)",
    nickel = "https://nickel-lang.org/",
    source = "stdlib/ByteDetector.affine (canonical pattern definitions)"
  }
}
