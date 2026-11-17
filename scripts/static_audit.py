#!/usr/bin/env python3
"""Non-network static gate for a clean GitHub release. Does NOT validate R runtime."""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
RESEARCH_EXT = {'.csv', '.tsv', '.xlsx', '.xls', '.rds', '.RData', '.fasta', '.fa', '.fastq',
                '.fq', '.bam', '.vcf', '.gpkg', '.sqlite', '.db', '.parquet', '.pdf', '.png', '.tiff'}
RESTRICTED_DIRS = {'data', 'results', 'private', 'raw', 'cache'}
CODE_EXT = {'.R', '.Rmd', '.py'}
SUSPICIOUS = [
    ("home_directory", re.compile(r'/(?:Users|home)/[^/\s\"\']+/')),
    ("literal_credentials", re.compile(r'(?:api[_-]?key|access[_-]?token|private[_-]?key|password|secret)\s*(?:<-|=|:)\s*[\"\'][^\"\']{6,}[\"\']', re.I)),
    ("credential_token", re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9_]{16,}|AKIA[A-Z0-9]{16}|sk-[A-Za-z0-9_-]{20,})\b')),
    ("sensitive_matching_code", re.compile(r'(?:by\s*=\s*[\"\']name[\"\']|by\s*=\s*[\"\']姓名[\"\'])')),
]


def r_delimiter_check(text):
    """Check matched (), {}, [] while ignoring comments/strings/backticks.

    This is intentionally NOT a replacement for parse(file=...) in R.
    """
    stack = []
    closing = {')': '(', '}': '{', ']': '['}
    quote = None
    comment = False
    escaped = False
    for lineno, line in enumerate(text.splitlines(), 1):
        comment = False
        for ch in line + '\n':
            if comment:
                continue
            if quote:
                if escaped:
                    escaped = False
                elif ch == '\\':
                    escaped = True
                elif ch == quote:
                    quote = None
                continue
            if ch == '#':
                comment = True
            elif ch in ('"', "'", '`'):
                quote = ch
            elif ch in '({[':
                stack.append((ch, lineno))
            elif ch in ')}]':
                if not stack or stack[-1][0] != closing[ch]:
                    return f"line {lineno}: unmatched {ch}"
                stack.pop()
    if quote:
        return "unclosed string/backtick"
    if stack:
        return f"unclosed {stack[-1][0]} from line {stack[-1][1]}"
    return None


def main():
    problems = []
    files = []
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '.git' in p.parts:
            continue
        rel = p.relative_to(ROOT)
        files.append(rel)
        if rel.parts[0] in RESTRICTED_DIRS:
            problems.append(f"Research/private directory must be empty: {rel}")
        if p.suffix.lower() in {e.lower() for e in RESEARCH_EXT}:
            problems.append(f"Disallowed data/output/binary extension: {rel}")
        if p.suffix not in CODE_EXT and p.suffix not in {'.md', '.gitignore'} and p.name != '.gitignore':
            problems.append(f"Unreviewed file type: {rel}")
        if p.suffix in CODE_EXT or p.suffix == '.md':
            content = p.read_text(encoding='utf8')
            for kind, pattern in SUSPICIOUS:
                if pattern.search(content):
                    problems.append(f"{kind}: {rel}")
            if p.suffix == '.R':
                issue = r_delimiter_check(content)
                if issue:
                    problems.append(f"R delimiter check: {rel} {issue}")
    required = {'.gitignore', 'README.md', 'docs/CODE_AUDIT.md', 'docs/PRIVACY_AND_RELEASE.md',
                'R/network_core.R', 'scripts/01_reconstruct.R', 'tests/test_core.R'}
    missing = sorted(required - {str(x) for x in files})
    for m in missing:
        problems.append('Missing essential file: ' + m)
    print(f'Checked {len(files)} files in cleaned project (no data read, no network access).')
    if problems:
        print('FAILED STATIC AUDIT:')
        for q in problems: print(' -', q)
        return 1
    print('PASS: no files in private directories, no banned research-file formats,')
    print('      no matched credential/name/home-directory patterns, balanced R delimiters.')
    print('NOT TESTED: R parser, R package versions, 10k network reconstruction,')
    print('            data-dependent numerical concordance, institutional clearance.')
    return 0

if __name__ == '__main__':
    sys.exit(main())
