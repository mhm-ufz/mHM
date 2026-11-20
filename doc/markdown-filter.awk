# Translate Markdown math and preserve percent signs in inline code for Doxygen.
# Keep one output line per input line so Doxygen can report source locations.

# Use matching backtick runs to find the end of a single-line code span.
function code_span_end(line, width,    offset) {
    offset = 0
    while (match(line, /`+/)) {
        offset += RSTART + RLENGTH - 1
        if (RLENGTH == width)
            return offset
        line = substr(line, RSTART + RLENGTH)
    }
    return 0
}

function print_inline(line,    token, width, end) {
    while (match(line, /`+|\\f\$([^$\\]|\\.)*\\f\$|\\.|\$\$+|\$(\\.|[^$\\[:space:]])([^$\\]|\\.)*\$/)) {
        token = substr(line, RSTART, RLENGTH)
        printf "%s", substr(line, 1, RSTART - 1)
        line = substr(line, RSTART + RLENGTH)

        if (token ~ /^`/) {
            width = length(token)
            end = code_span_end(line, width)
            if (end) {
                token = token substr(line, 1, end)
                line = substr(line, end + 1)
                # Doxygen consumes the percent sign immediately before a name.
                gsub(/%[[:alpha:]_]/, "%&", token)
            }
        } else if (token ~ /^\$[^$]/ && token !~ /[[:space:]]\$$/ && line !~ /^[[:digit:]]/) {
            token = "\\f" substr(token, 1, length(token) - 1) "\\f$"
        }
        printf "%s", token
    }
    print line
}

function print_display(    environment, closing, i) {
    environment = display[2]
    if (environment ~ /^[[:blank:]]*\\begin\{(equation|align)\*?\}[[:blank:]]*$/) {
        sub(/^[[:blank:]]*\\begin\{/, "", environment)
        sub(/\}[[:blank:]]*$/, "", environment)
    } else {
        environment = ""
    }
    closing = display[n_display - 1]
    sub(/^[[:blank:]]*/, "", closing)
    sub(/[[:blank:]]*$/, "", closing)
    if (environment && closing == "\\end{" environment "}") {
        print ""
        print "\\f{" environment "}{"
        for (i = 3; i < n_display - 1; i++)
            print display[i]
        print "\\f}"
        print ""
    } else {
        print "\\f["
        for (i = 2; i < n_display; i++)
            print display[i]
        print "\\f]"
    }
}

{
    # A closing fence must match the kind and be at least as long as the opener.
    fence = ""
    if (!n_display && match($0, /^ ? ? ?(```+|~~~+)/)) {
        fence = substr($0, 1, RLENGTH)
        sub(/^ +/, "", fence)
        after_fence = substr($0, RLENGTH + 1)
    }
    if (code_fence || (fence && (fence ~ /^~/ || after_fence !~ /`/))) {
        if (!code_fence)
            code_fence = fence
        else if (substr(fence, 1, 1) == substr(code_fence, 1, 1) &&
                 length(fence) >= length(code_fence) && after_fence ~ /^[[:blank:]]*$/)
            code_fence = ""
        print
    } else if ($0 ~ /^ ? ? ?\$\$[[:blank:]]*$/) {
        display[++n_display] = $0
        if (n_display > 1) {
            print_display()
            n_display = 0
        }
    } else if (n_display) {
        display[++n_display] = $0
    } else {
        print_inline($0)
    }
}

# Leave an unmatched display delimiter and its trailing text as written.
END {
    for (i = 1; i <= n_display; i++)
        print display[i]
}
