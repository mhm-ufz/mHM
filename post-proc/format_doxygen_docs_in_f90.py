#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : format_doxygen_docs_in_f90
Project Name: mpr_extract
Description : insert your description here, if applicable
Author      : ottor
Created     : 16.05.18 15:23
"""

# IMPORTS
import pathlib
import re
import pytest
from datetime import datetime

# GLOBAL VARIABLES
POSSIBLE_SUFFIXES = ['.f90', '.F']
ROUTINE_FLAG = ['function', 'subroutine']
MODIFIED_FLAG = ['Modified', 'Modifications']
IMPLICIT_NONE_FLAG = 'implicit none'
USE_FLAG = r'^[\s]*use'
ARG_DESCRIPTION_SPLITTER = '::'
MAX_LINE_WIDTH = 100

# FUNCTIONS
def get_all_subfiles(path, relation=None):
    if path.is_file() and path.suffix in POSSIBLE_SUFFIXES:
        yield path.relative_to(relation or path)
    else:
        for sub_path in path.iterdir():
            yield from get_all_subfiles(sub_path, path)

def check_for_line_to_ignore(line):
    return not line.strip() or re.match('!\s+-+\s*$', line.strip())

# CLASSES
class Doc(object):
    FILE_FLAGS = ['brief', 'details', 'authors', 'author', 'date', 'note', 'return']
    BLOCK_DICT = {'brief': 'PURPOSE', 'authors': 'HISTORY', 'note':'RESTRICTIONS',
                  'return': 'RETURN'}
    DOC_COMMENT = '!>'
    COMMENT = '!'
    FLAG_PREFIX = '\\'

    def __init__(self):
        self._name = None
        self.param_dict = {}
        self._brief = None
        self._details = None
        self._return = None
        self._author = None
        self._authors = None
        self._note = None
        self._example = None
        self._literature = None
        self._date = None

        self._in_modified = False
        self._modified = []
        self.cur_field = None

    def parse_doc(self, line):
        if not check_for_line_to_ignore(line):
            if any([re.search(pattern, line) for pattern in MODIFIED_FLAG]):
                self._in_modified = True

            if self._in_modified:
                self.parse_modifications(line)
            elif line.strip().startswith(self.DOC_COMMENT):
                # check for DOC_COMMENT
                parsed = None
                for pattern in self.FILE_FLAGS:
                    matched = re.match(self.DOC_COMMENT +
                             '\s*' +
                             self.FLAG_PREFIX * 2 +
                             pattern, line.strip(), re.IGNORECASE)
                    if matched:
                        self.cur_field = pattern
                        parsed = line.strip()[matched.end():].strip()
                        self.set_file_doc_attr(parsed)
                if parsed is None:
                    self.set_file_doc_attr(line.lstrip(' ' + self.DOC_COMMENT))

    def doc_to_str(self, fname=None, indent=0, module=False):

        def format_item(indent, comment, block_name, value):
            value = value.replace('\n', '\n{indent}{comment:9}'.format(indent=' ' * indent,
                                                                       comment=comment))
            template = '{indent}{comment:5}{block_name}\n{indent}{comment:9}{value}\n\n'
            if block_name is None:
                template = '{indent}{comment:9}{value}\n\n'
            return template.format(indent=' ' * indent,
                                   comment=comment,
                                   block_name=block_name,
                                   value=value)

        value = getattr(self, '_author')
        if value is not None:
            setattr(self, '_author', None)
            setattr(self, '_authors', value)

        out_string = ''

        if module:
            out_string += format_item(indent=indent, comment=self.COMMENT,
                                      block_name='NAME', value=self._name)
        for item in self.FILE_FLAGS:
            value = getattr(self, '_{}'.format(item))
            if item == 'file' and fname is not None:
                value = fname
            if value is None:
                if item == 'authors':
                    value = self.default_author
                elif item == 'date':
                    value = '{:%b %Y}'.format(datetime.today())
                elif item in ['brief', 'details']:
                    value = 'TODO: add description'
            if value is not None:
                value = '{} {}'.format(self.FLAG_PREFIX + item, value)
                out_string += format_item(indent=indent, comment=self.DOC_COMMENT,
                                          block_name=self.BLOCK_DICT.get(item), value=value)

        out_string += '{}{} Modifications:\n'.format(' ' * indent, self.COMMENT)
        if self._modified:
            max_editor_length = max((mod.editor_length for mod in self._modified))
            for mod in self._modified:
                mod.note = '- ' + mod.note.lstrip('- ').replace('- ', '\n{}{}{}- '.format(' ' * (indent),
                                                                                          self.COMMENT,
                                                                                          ' ' * (11 + max_editor_length)))
                out_string += '{}{} {:{width}} {} {}\n'.format(' ' * indent, self.COMMENT,
                                                       mod.editor, mod.date, mod.note, width=max_editor_length)

        return out_string + '\n'

    def parse_modifications(self, line):
        split = re.split('([A-Z][a-z]{2}\s+20\d\d)',line)
        if len(split) == 3:
            # modified is in line
            for pattern in MODIFIED_FLAG:
                matched = re.search(pattern, split[0])
                if matched:
                    split[0] = split[0][matched.end():]
            editor = split[0].strip('! ')
            date = split[1].strip()
            note = split[2].strip()
            self._modified.append(Modification(editor, date, note))
        else:
            self._modified[-1].note += ' {}'.format(line.strip('! \n'))

    def set_file_doc_attr(self, value):
        if self.cur_field is None:
            raise Exception
        cur_name = '_{}'.format(self.cur_field)

        cur_attr = getattr(self, cur_name)
        value = value.strip().replace('\\n', '')
        if cur_attr is None:
            # get the info and append to string
            setattr(self, cur_name, value)
        else:
            setattr(self, cur_name, cur_attr + '\n' + value)


class FortranFile(Doc):
    FILE_FLAGS = ['file', 'brief', 'details', 'author', 'authors', 'date', 'version', 'copyright']
    BLOCK_DICT = {}

    MODULE_FLAG = ['module', 'program']

    def __init__(self, default_author):
        # parameters
        super().__init__()
        self._file = None
        self._version = None
        self._copyright = None
        self.lines = []
        self.comments = None
        self.modifications = []
        self._cur_doc = None

        # flags
        self._cur_routine = None
        self._in_module = False
        self._in_modified = False
        self._in_doc = False
        self._in_routine = False
        self._in_routine_args = False
        self._is_line_continued = False
        self._in_code = False

        self.default_author = default_author

    def read(self, fname):
        print('reading', fname)
        with open(fname, 'r') as f_in:
            for line in f_in:
                self._check_status(line)

                if not self._in_module:
                    self._parse_file_doc(line)
                elif self._in_doc:
                    self._parse_doc(line)
                elif self._in_routine:
                    self._parse_routine(line)
                elif self._in_code:
                    if self._cur_routine is not None:
                        self.lines.append(self._cur_routine)
                        if self._cur_routine.comment_cache:
                            self.lines.extend(self._cur_routine.comment_cache)
                        self._cur_routine = None
                    self.lines.append(line)
                # elif line.strip():
                else:
                    self.lines.append(line)

    def write(self, fname):
        print('writing', fname)
        with open(fname, 'w') as f_in:
            f_in.write(self.doc_to_str(fname.name))
            for line in self.lines:
                if isinstance(line, Routine):
                    f_in.write(line.to_str())
                else:
                    f_in.write(line)

    def _check_status(self, line):
        if check_for_line_to_ignore(line):
            pass
        elif not self._in_module:
            self._check_for_file_doc(line)
        elif self._in_doc:
            self._check_for_routine(line)
        elif self._in_routine_args:
            pass
        elif self._in_routine:
            self._check_for_code(line)
        elif self._in_code:
            self._check_for_code_end(line)
        else:
            self._check_for_doc(line)
            self._check_for_routine(line)

    def _check_for_file_doc(self, line):
        # check for MODULE_FLAG???
        if any([re.search(pattern, line.split(self.COMMENT)[0], re.IGNORECASE) for pattern in self.MODULE_FLAG]):
            self._in_module = True
            self._in_modified = False
            self._cur_filedoc = None

    def _check_for_routine(self, line):
        # check for ROUTINE_FLAG
        if any([re.search(pattern, line.split(self.COMMENT)[0], re.IGNORECASE) for pattern in ROUTINE_FLAG]):
            self._in_doc = False
            self._in_modified = False
            self._in_routine = True
            self._in_routine_args = True

    def _check_for_line_continuation(self, line):
        stripped = line.split(self.COMMENT)[0].strip()
        if stripped.endswith('&'):
            self._is_line_continued = True
        elif stripped:
            self._is_line_continued = False

    def _check_for_doc(self, line):
        # check for DOC_COMMENT
        if re.match(self.DOC_COMMENT, line.strip()) or re.match('!\s+NAME', line.strip()):
            self._in_doc = True

    def _check_for_code(self, line):
        # check for neither use nor "implicit none" nor :: and line
        clean_line = line.split(self.COMMENT)[0].strip()
        if not (ARG_DESCRIPTION_SPLITTER in clean_line or
                re.match(USE_FLAG, clean_line, re.IGNORECASE) or
                re.match(IMPLICIT_NONE_FLAG, clean_line, re.IGNORECASE)) \
                and not self._is_line_continued and clean_line:
            self._in_routine = False
            self._in_code = True

    def _check_for_code_end(self, line):
        # check for end ROUTINE_FLAG
        if any([re.match('end\s*' + pattern, line.split(self.COMMENT)[0], re.IGNORECASE) for pattern in ROUTINE_FLAG]):
            self._in_code = False

    def _parse_doc(self, line):
        if self._cur_doc is None:
            self._cur_doc = Doc()
        # check for DOC_COMMENT
        self._cur_doc.parse_doc(line)

    def _parse_file_doc(self, line):
        # update self._cur_doc
        self.parse_doc(line)

    def _parse_routine(self, line):
        if self._cur_routine is None:
            self._cur_routine = Routine(doc=self._cur_doc)
            # update self._cur_doc
            self._cur_doc = None
        self._cur_routine.parse_line(line,
                                     is_line_continued=self._is_line_continued,
                                     in_routine_args=self._in_routine_args,
                                     )
        self._check_for_line_continuation(line)
        if not self._is_line_continued and self._in_routine_args:
            self._in_routine_args = False


class Routine(object):
    ONLY_STR = 'only'
    COMMENT = '!'
    def __init__(self, doc=None):
        self.doc = doc
        self.args = []
        self.args_doc = []
        self.args_attrs = []
        self.implicit_none = False
        self.used = []
        self.name = None
        self.type = None
        self.type_attrs = None
        self.type_addon = None
        self.comment_cache = []

    def parse_line(self, line, is_line_continued=False, in_routine_args=False):
        # check for DOC_COMMENT
        if self.name is None:
            for pattern in ROUTINE_FLAG:
                matched = re.match(pattern, line.strip(), re.IGNORECASE)
                if not matched:
                    continue
                # get the type
                self.type = pattern
                # get the stuff like (elemental pure...)
                self.type_attrs = line.strip()[:matched.start()] or None
                # trim line
                line = line.strip()[matched.end():]
                split = line.split('(')
                self.name = split[0].strip()
                line = split[1]
                in_routine_args = True
            if self.name is None:
                print(line)
                raise Exception
        if in_routine_args:
            # remove comments
            split = line.split('!')
            if len(split) == 2:
                comment = split[1]
            else:
                comment = ''
            split = split[0].split(')')
            if len(split) > 1 and split[1].rstrip('& \n'):
                self.type_addon = split[1].rstrip('& \n') + ')'
            cur_args = [arg.strip(' &') for arg in split[0].split(',') if arg.strip(' &')]
            if len(cur_args) == 1 and comment:
                self.args_doc += [comment]
            else:
                self.args_doc += [comment] * len(cur_args)
            self.args += cur_args
        else:
            self._parse_body(line, is_line_continued)

    def _parse_body(self, line, is_line_continued):
        if not line.strip().startswith(self.COMMENT):
            self.comment_cache = []
        matched = re.match(IMPLICIT_NONE_FLAG, line.strip(), re.IGNORECASE)
        if matched:
            self.implicit_none = True
        elif re.match(USE_FLAG, line.strip(), re.IGNORECASE):
            line = line.strip()[re.match(USE_FLAG, line.strip(), re.IGNORECASE).end():]
            matched = re.search(self.ONLY_STR, line, re.IGNORECASE)
            if not matched:
                raise Exception('use without "{}" discovered'.format(self.ONLY_STR))
            split = re.split(self.ONLY_STR, line, flags=re.IGNORECASE)
            split[0] = split[0].strip(' ,').lower()
            for item in split[1].split(','):
                value = item.split(self.COMMENT)[0].strip(':&\n ')
                if value:
                    self.used.append((split[0], value))
        elif ARG_DESCRIPTION_SPLITTER in line:
            split = line.split(ARG_DESCRIPTION_SPLITTER)
            name = split[1].split('!')
            self.args_attrs.append(ArgumentAttributes(split[0], *name))
        elif line.strip().startswith('!'):
            if line.count(self.COMMENT) >= 2:
                self.args_attrs[-1].append_str(line.lstrip(self.COMMENT + ' '))
            #elif not check_for_line_to_ignore(line):
            else:
                self.comment_cache.append(line)
                # self.args_attrs.append(line.lstrip('! '))
                # raise Exception
        elif is_line_continued:
            module = self.used[-1][0]
            for item in line.split(','):
                self.used.append((module, item.strip(':&\n ')))


    def to_str(self):
        routine_str = ''
        indent=2
        # TODO: alter the doc depending on args
        self.update_doc()
        routine_str += self.doc.doc_to_str(indent=indent, module=True)
        routine_str += self._args_to_str(indent=indent)
        routine_str += self._imports_to_str(indent=indent+2)
        routine_str += self._implicit_to_str(indent=indent+2)
        routine_str += self._variables_to_str(indent=indent+2)
        return routine_str

    def _args_to_str(self, indent=0):
        args_str = '{}{}{}('.format(' ' * indent,
                                 self.type_attrs or '',
                                   self.type,
                                   )
        arg_length = len(args_str)
        total_length = 0
        for arg in self.args:
            if arg_length + total_length + len(arg) + 2 > MAX_LINE_WIDTH:
                args_str += '&\n{}'.format(' ' * arg_length)
                total_length = 0
            args_str += '{}, '.format(arg)
            total_length += len(arg) + 2

        return '{}){}\n\n'.format(args_str[:-2], self.type_addon or '')

    def _imports_to_str(self, indent=0):
        def format_use(use_item):
            new_string = '{}use {}, {} {} {}'.format(' ' * indent,
                                             use_item[0],
                                             self.ONLY_STR,
                                             ':',
                                             use_item[1])
            return new_string, new_string.find(':') + 2
        if self.used:
            used = sorted(self.used)
            out_string, arg_length = format_use(used[0])
            total_length = 0
            for i_item, use_item in enumerate(used[1:], 1):
                if use_item[0] == used[i_item-1][0]:
                    if arg_length + total_length + len(use_item[1]) + 2 > MAX_LINE_WIDTH:
                        out_string += ', &\n{}{}'.format(' ' * arg_length, use_item[1])
                        total_length = len(use_item[1])
                    else:
                        out_string += ', {}'.format(use_item[1])
                        total_length += len(use_item[1]) + 2
                else:
                    out_string_addon, arg_length = format_use(use_item)
                    out_string += '\n' + out_string_addon
            return '{}\n\n'.format(out_string)
        else:
            return ''

    def _variables_to_str(self, indent):
        out_string = ''
        for attr in self.args_attrs:
            comment = ''
            if attr.doc:
                comment = ' {} {}'.format(self.COMMENT, attr.doc)
            out_string += '{}{} {} {}{}\n'.format(' ' * indent,
                                             ', '.join(attr.attrs),
                                             ARG_DESCRIPTION_SPLITTER,
                                             attr.name,
                                             comment)
        return out_string + '\n'

    def _implicit_to_str(self, indent):
        if not self.implicit_none:
            raise Exception('implicit none is not set in routine' + self.name)
        return '{}{}\n\n'.format(' ' * indent, IMPLICIT_NONE_FLAG)

    def update_doc(self):
        self.doc._name = self.name
        for arg, arg_doc in zip(self.args, self.args_doc):
            arg_attrs = self._get_args_attr(arg)
            # set the doc it it was only in set in argument list
            if not arg_attrs.doc and arg_doc:
                arg_attrs.doc = arg_doc
            if arg_attrs.intent_optional_key:
                self.doc.param_dict.setdefault(arg_attrs.intent_optional_key, []).append(arg_attrs.doxygen_string)


    def _get_args_attr(self, name):
        for arg in self.args_attrs:
            if arg.name == name:
                return arg
        raise KeyError('variable ' + name + ' was not found in the AttributeTable')


class ArgumentAttributes(object):
    def __init__(self, attrs, name, docstring=''):
        split = name.strip().split('(')
        self.name = split[0].strip()
        self.attrs = []
        self.intent_attr = ''
        self.optional_attr = ''
        self.dimension_attr = ''

        for i_attr, attr in enumerate(attrs.split(',')):
            self.attrs.append(attr.strip().lower())
            if i_attr == 0:
                self.type_attr = attr.strip(' ')
            if attr.startswith('intent'):
                self.intent_attr = attr.split('(')[1].strip(') ')
            if attr.startswith('optional'):
                self.optional_attr = attr.strip(' ')
            if attr.startswith('dimension'):
                self.dimension_attr = attr.strip(' ')

        if len(split) == 2:
            attr_str = 'dimension(' + split[1]
            self.attrs.append(attr_str)
            self.dimension_attr = attr_str
        self.doc = docstring.strip() or ''

    def append_str(self, to_append):
        self.doc += to_append.strip()

    @property
    def intent_optional_key(self):
        keys = []
        keys += [self.intent_attr]
        keys += [self.optional_attr]
        return ''.join(keys)

    @property
    def doxygen_string(self):
        param_args = [arg for arg in [self.type_attr, self.dimension_attr] if arg]
        return ('\\param[{}] "{} {} {}"'.format(self.intent_attr,
                                                ', '.join(param_args),
                                                ARG_DESCRIPTION_SPLITTER,
                                                self.name),
                self.doc)
        """e.g. \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with"""

class Modification(object):
    def __init__(self, editor, date, note):
        self.editor = editor.strip(', ')
        self.date = date
        self.note = note

    @property
    def editor_length(self):
        return len(self.editor)

class TestFortranFile(object):
    @pytest.fixture
    def fortran_file(self):
        return FortranFile()

    def test_check_for_file_doc(self, fortran_file):
        lines = ['module example(&', 'MODULE EXAMPLE(&', 'module example(a, b, c)']
        for line in lines:
            fortran_file._in_module = False
            fortran_file._check_for_file_doc(line)
            assert fortran_file._in_module

    def test_check_for_routine(self, fortran_file):
        lines = ['function newDimAlias(aliases) result(newDimAlias)',
                 'SUBROUTINE mpr_eval(parameterset)',
                 'function majority_statistics(nClass, &   ! number of classes',
                 '  elemental pure subroutine temporal_disagg_forcing(isday, &']
        for line in lines:
            fortran_file._in_routine = False
            fortran_file._check_for_routine(line)
            assert fortran_file._in_routine

    def test_check_for_doc(self, fortran_file):
        lines = ['  !>        \param[out] "real(dp), dimension(:) :: aet"                 actual ET [mm/s]',
                 '!>        pervious areas is calculated as (omit \f$t\f$)',
                 '  !>        \param[in] "integer(i4),           :: processCase"']
        for line in lines:
            fortran_file._in_doc = False
            fortran_file._check_for_doc(line)
            assert fortran_file._in_doc

    def test_check_for_code(self, fortran_file):
        lines = ['    runoff_sealed = 0.0_dp',
                 '    end do ! hh']
        for line in lines:
            fortran_file._in_code = False
            fortran_file._check_for_code(line)
            assert fortran_file._in_code

    def test_check_for_code_end(self, fortran_file):
        lines = ['END   FUNCTION',
                 '    end subroutine dummy']
        for line in lines:
            fortran_file._in_code = True
            fortran_file._check_for_code_end(line)
            assert not fortran_file._in_code


# SCRIPT
if __name__ == '__main__':
    default_author = 'Robert Schweppe'
    arg_path_in = '/Users/ottor/ownCloud/Home/local_libs/fortran/mpr_extract/src/mHM/'
    arg_path_in = '/Users/ottor/ownCloud/Home/local_libs/fortran/mpr_extract/src/mHM/mo_mhm.f90'
    arg_path_out = '/Users/ottor/temp/mpr_extract_test/'
    path_in = pathlib.Path(arg_path_in)

    if path_in.is_file():
        path_list = [path_in.name]
        path_in = path_in.parent
    else:
        path_list = get_all_subfiles(path_in)

    for path in path_list:
        filereader = FortranFile(default_author)
        filereader.read(pathlib.Path(path_in, path))
        filereader.write(pathlib.Path(arg_path_out, path))