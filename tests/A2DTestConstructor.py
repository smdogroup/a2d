import datetime
import inspect
import itertools
import os
import re
import warnings
from typing import List, Dict, Callable, Any, Set, Tuple, Union, FrozenSet
import numpy

__package_name__ = __package__ if __package__ else "A2DTestConstructor"
generator_return_type = numpy.ndarray
var_generator_type = Callable[[], generator_return_type]


class VarType:
    """\
VarType: A representation of a constant or non-differentiated variable type.  This object contains all the \
information necessary to implement use of the type in a testing script.

alias: the name to be defined and used by the C++ program for the type (usually ending in "_t").

shortname: an abbreviated name for the type, used when describing the input types in the "TEST_F" parent classes.

definition: the specific C or C++ type used to define the type as alias.  For example, "A2D::Vec<T, 3>" or "double".

test: the templated test function used by C++ to test if two values of this type are equal.  For example, \
"expect_val_eq" for scalar types or "expect_vec_eq<3>" for length 3 vector types.

initialize_assignment_format: a format string used in conjunction with the "generator" to define a random instance \
of the type.  This should usually have the form "const <alias> {var_name} = {0};" for scalar types or \
"const <type> {var_name}[N] = {{{0}, {1}, {2}, ..., {N}}};" for iterable types, where <type> represents the base type \
and N is the iterable length.

rank: the integer rank of the object.  The rank is equal to the number of indices used to exactly specify a single \
element of the object.  For example, scalars are rank 0, vectors are rank 1, matrices are rank 2, etc...

generator: a zero-argument function returning a numpy array of values.  The function must return a numpy array \
in order to be compatible with the complex step method implementations.  The values returned by the "generator" \
function should be sufficient to initialize a random instance of the variable type when formatted as a variadic \
argument to the "initialize_assignment_format" string using the str.format method.  For example, \
self.initialize_format.format(*self.generator(), **kwargs)
"""

    __instances_by_alias__: Dict[str, "VarType"] = {}

    alias: str
    shortname: str
    definition: str
    test: str
    initialize_assignment_format: str
    rank: int
    generator: var_generator_type

    def __init__(self,
                 alias: str,
                 shortname: str,
                 definition: str,
                 test: str,
                 initialize_assignment_format: str,
                 rank: int,
                 generator: var_generator_type,
                 ):
        self.alias = alias
        self.shortname = shortname
        self.definition = definition
        self.test = test
        self.initialize_assignment_format = initialize_assignment_format
        self.rank = rank
        self.generator = generator

        self.alias_pattern = re.compile(fr'(\A|[., <>()\[\]]){self.alias}([., <>()\[\]]|\Z)')
        pass

    def __new__(cls,
                alias: str,
                shortname: str,
                definition: str,
                test: str,
                initialize_assignment_format: str,
                rank: int,
                generator: var_generator_type,
                ):
        if alias in cls.__instances_by_alias__:
            equivalence = [getattr(cls.__instances_by_alias__[alias], key) == locals()[key]
                           for key in cls.__comparison_keys__()]
            if all(equivalence):
                # check generator
                try:
                    # noinspection PyUnresolvedReferences
                    generator_check = (cls.__instances_by_alias__[alias].generator.__code__.co_code ==
                                       generator.__code__.co_code)
                except AttributeError as e:
                    w = UserWarning(f'\nThe check for agreement between generator objects for the VarType'
                                    f'\nobject "{alias}" did not succeed.  There is enough agreement'
                                    f'\nbetween other attributes to assume the definitions are equivalent.'
                                    f'\nNo new object was created, instead the previously defined "{alias}"'
                                    f'\nobject has been returned.')
                    warnings.warn(w)
                    generator_check = True
                if generator_check:
                    return cls.__instances_by_alias__[alias]
            raise ValueError(f'Two different VarType objects with the same alias are not allowed.  There is already\n'
                             f'a VarType object with the alias "{alias}" with different attributes than were defined.')
        if alias in ADVarType.__instances_by_alias__:
            raise ValueError(f'A VarType object cannot have the same alias as a ADVarType object.')
        if alias in A2DVarType.__instances_by_alias__:
            raise ValueError(f'A VarType object cannot have the same alias as a A2DVarType object.')

        obj = object.__new__(cls)
        cls.__instances_by_alias__[alias] = obj
        return obj

    def __eq__(self, other):
        if isinstance(other, VarType):
            equivalence = all(getattr(self, key) == getattr(other, key)
                              for key in self.__comparison_keys__())
            return equivalence
        return False

    def __hash__(self):
        return sum(map(hash, (
            getattr(self, key)
            for key in self.__comparison_keys__().union(['generator'])
        )))

    @classmethod
    def __comparison_keys__(cls):
        return set(cls.__init__.__annotations__.keys()).difference(['generator'])

    pass


class ADVarType:
    """\
ADVarType: A representation of an automatically differentiated variable type.  This object contains all the \
information necessary to implement use of the type in a testing script.

alias: the name to be defined and used by the C++ program for the type (usually ending in "_t").

shortname: an abbreviated name for the type, used when describing the input types in the "TEST_F" parent classes.

definition: the specific AD type used to define the type as alias.  For example, "A2D::ADVec<Vec_t>" or \
"A2D::ADScalar<T>".

get_value: the method call operation used to obtain the base value for the type.  Usually ".value" or ".value()".

get_bvalue: the method call operation used to obtain the seed or bar value for the type.  Usually ".bvalue" or \
".bvalue()".

parent: the VarType object this ADVarType object inherits from.  Both the value and bvalue for objects of this type \
must be of the parent type, and initialization of this type must take the form "<alias> <name>(<p>, <pb>);" where \
p and pb are instances of the parent object.

initialization_defaults: (optional) values for initialization of a new instance of this type, for use when declaration \
of the "value" and "bvalue" parent objects is not necessary for the AD class to be instantiated.  MUST be defined when \
such is the case, but should be left missing when not.  This is usually necessary when the parent type is not passed \
to the constructor by reference.
"""
    __instances_by_alias__: Dict[str, "ADVarType"] = {}

    alias: str
    shortname: str
    definition: str
    get_value: str
    get_bvalue: str
    parent: VarType
    initialization_defaults: Union[Tuple[Any, Any], Tuple[()]]

    # NOTE: this REQUIRES that this type be constructed as:
    #   <alias> <name>(<parent_obj>, <parent_b_obj>);
    # EX:
    #   ADVec_t x_ad(x, xb);
    # Where x and xb are Vec_t types

    def __init__(self,
                 alias: str,
                 shortname: str,
                 definition: str,
                 get_value: str,
                 get_bvalue: str,
                 parent: VarType,
                 initialization_defaults: Union[Tuple[Any, Any], Tuple[()]] = (),
                 ):
        self.alias = alias
        self.shortname = shortname
        self.get_value = get_value
        self.get_bvalue = get_bvalue
        self.parent = parent
        self.initialization_defaults = initialization_defaults

        # keep from overwriting the (dynamic) definition property of child class(es)
        if not hasattr(self.__class__, 'definition'):
            self.definition = definition
        self.alias_pattern = re.compile(fr'(\A|[., <>()\[\]]){self.alias}([., <>()\[\]]|\Z)')
        pass

    def __new__(cls,
                alias: str,
                shortname: str,
                definition: str,
                get_value: str,
                get_bvalue: str,
                parent: VarType,
                initialization_defaults: Union[Tuple[Any, Any], Tuple[()]] = (),
                ):
        if alias in cls.__instances_by_alias__:
            equivalence = all(getattr(cls.__instances_by_alias__[alias], key) == locals()[key]
                              for key in cls.__comparison_keys__())
            if equivalence:
                return cls.__instances_by_alias__[alias]
            raise ValueError(f'Two different ADVarType objects with the same alias are not allowed.  There is already\n'
                             f'an ADVarType object with the alias "{alias}" with different attributes than were given.')
        if alias in VarType.__instances_by_alias__:
            raise ValueError(f'An ADVarType object cannot have the same alias as a VarType object.')
        if alias in A2DVarType.__instances_by_alias__:
            raise ValueError(f'An ADVarType object cannot have the same alias as an A2DVarType object.')

        obj = object.__new__(cls)
        cls.__instances_by_alias__[alias] = obj
        return obj

    @property
    def test(self):
        return self.parent.test

    @property
    def initialize_assignment_format(self):
        return self.parent.initialize_assignment_format

    @property
    def rank(self):
        return self.parent.rank

    def __eq__(self, other):
        if isinstance(other, ADVarType):
            equivalence = all(getattr(self, key) == getattr(other, key)
                              for key in self.__comparison_keys__())
            return equivalence
        return False

    def __hash__(self):
        return sum(map(hash, (
            getattr(self, key)
            for key in self.__comparison_keys__()
        )))

    @classmethod
    def __comparison_keys__(cls):
        return cls.__init__.__annotations__.keys()

    @property
    def requires_parent_declaration(self):
        return self.initialization_defaults == ()

    pass


class A2DVarType(ADVarType):
    """\
A2DVarType: A representation of a doubly automatically differentiated variable type.  This object contains all the \
information necessary to implement use of the type in a testing script.

alias: the name to be defined and used by the C++ program for the type (usually ending in "_t").

shortname: an abbreviated name for the type, used when describing the input types in the "TEST_F" parent classes.

definition: the specific A2D type used to define the type as alias.  The "N" template parameter must be included as \
a string format option ("{N}") in the definition.  For example, "A2D::A2DVec<{N}, Vec_t>" or "A2D::A2DScalar<{N}, T>".

get_value: the method call operation used to obtain the base value for the type.  Usually ".value" or ".value()".

get_bvalue: the method call operation used to obtain the seed or bar value for the type.  Usually ".bvalue" or \
".bvalue()".

get_pvalue: the method call operation formatted to obtain the seed or dot value for the type from an array of such \
values.  This must be a format string with variadic inputs.  Usually ".pvalue[{0}]" or ".pvalue({0})".

get_hvalue: the method call operation formatted to obtain the seed or hat value for the type from an array of such \
values.  This must be a format string with variadic inputs.  Usually ".hvalue[{0}]" or ".hvalue({0})".

parent: the VarType object this A2DVarType object inherits from.  The value, bvalue, pvalues, and hvalues for objects \
of this type must be of the parent type, and initialization of this type must take the form \
"<alias> <name>(<p>, <pb>);" where p and pb are instances of the parent object.

initialization_defaults: (optional) values for initialization of a new instance of this type, for use when declaration \
of the "value" and "bvalue" parent objects is not necessary for the A2D class to be instantiated.  MUST be defined \
when such is the case, but should be left missing when not.  This is usually necessary when the parent type is not \
passed to the constructor by reference.
"""
    __instances_by_alias__: Dict[str, "A2DVarType"] = {}

    alias: str
    shortname: str
    __definition__: str
    N: int
    get_value: str
    get_bvalue: str
    get_pvalue: str
    get_hvalue: str
    parent: VarType
    initialization_defaults: Union[Tuple[Any, Any], Tuple[()]]

    # NOTE: this REQUIRES that this type be constructed as:
    #   <alias> <name>(<parent_obj>, <parent_b_obj>);
    # EX:
    #   ADVec_t x_ad(x, xb);
    # Where x and xb are Vec_t types

    # noinspection PyPep8Naming
    def __init__(self,
                 alias: str,
                 shortname: str,
                 definition: str,
                 get_value: str,
                 get_bvalue: str,
                 get_pvalue: str,
                 get_hvalue: str,
                 parent: VarType,
                 N: int,
                 initialization_defaults: Union[Tuple[Any, Any], Tuple[()]] = (),
                 ):
        if "{N}" not in definition:
            raise ValueError(f'{{N}} must be included in the "definition" attribute of all A2DVarType objects,\n'
                             f'in order to dynamically define the value of "N".  "N" is either statically defined '
                             f'within\n"{definition}", or not defined at all.')
        self.get_pvalue = get_pvalue
        self.get_hvalue = get_hvalue
        self.N = N
        # noinspection PyTypeChecker
        ADVarType.__init__(self, alias, shortname, None, get_value, get_bvalue, parent, initialization_defaults)
        self.__definition__ = definition
        pass

    # noinspection PyPep8Naming
    def __new__(cls,
                alias: str,
                shortname: str,
                definition: str,
                get_value: str,
                get_bvalue: str,
                get_pvalue: str,
                get_hvalue: str,
                parent: VarType,
                N: int,
                initialization_defaults: Union[Tuple[Any, Any], Tuple[()]] = (),
                ):
        if alias in cls.__instances_by_alias__:
            equivalence = all(getattr(cls.__instances_by_alias__[alias], key) == locals()[key]
                              for key in cls.__comparison_keys__())
            if equivalence:
                return cls.__instances_by_alias__[alias]
            raise ValueError(f'Two different A2DVarType objects with the same alias are not allowed.  There is \n'
                             f'already an A2DVarType object with the alias "{alias}" with different attributes than '
                             f'were given.')
        if alias in VarType.__instances_by_alias__:
            raise ValueError(f'An A2DVarType object cannot have the same alias as a VarType object.')
        if alias in ADVarType.__instances_by_alias__:
            raise ValueError(f'An A2DVarType object cannot have the same alias as an ADVarType object.')

        obj = object.__new__(cls)
        cls.__instances_by_alias__[alias] = obj
        return obj

    @property
    def definition(self):
        return self.__definition__.format(N=self.N)

    pass


class VariableOverloadError(Exception):
    """Multiple types (VarTypes, ADVarTypes, and/or A2DVarTypes) assigned to same variable name."""
    pass


class NamingConventionConflict(Exception):
    """
    Collision of variable names for set naming convention.
    Using a variable name "<var_name>" reserves all names with form
    "<var_name>", "<var_name>b", "<var_name>p<#>", and "<var_name>h<#>".
    """
    pass


class VariableNameOverloadError(Exception):
    """Multiple variables with the same name."""
    pass


io_ad_type = List[Tuple[str, ADVarType]]  # [(<var_name>, <var_type>), ...]
io_a2d_type = List[Tuple[str, A2DVarType]]  # [(<var_name>, <var_type>), ...]
non_constant_inputs_type = FrozenSet[str]  # frozenset(<var_name>, ...)  <-These are the active inputs
non_constant_inputs_input_type = Set[str]  # {<var_name>, ...}  <-These are the active inputs
variants_type = List[non_constant_inputs_type]
variants_input_type = Union[List[non_constant_inputs_input_type],
                            Set[non_constant_inputs_type],
                            variants_type]
input_data_type = Dict[str, Tuple[generator_return_type, ADVarType]]  # {<var_name>: (<value>, <var_type>), ...}
operation_type = Callable[..., Tuple[numpy.ndarray, ...]]


class TestADFunction:
    """\
TestADFunction: A representation of an overloaded operation to be tested.

operation_name: the name of the overloaded operation within the A2D C++ code as a string using double quotations.  \
The operation should return an instance of the resulting A2D class, specifically an "AD" class; all outputs and \
modified objects must be passed by reference to the operation.

inputs: a list of variable-name, variable-type pairs that define the inputs to the overloaded operation in the MOST \
differentiable state, that is, when no input is considered constant.  The variable names should be strings, and the \
variable types should be ADVarType objects.

outputs: a list of variable-name, variable-type pairs that define the outputs to the overloaded operation in the MOST \
differentiable state.  For most if not all circumstances, all outputs should be differentiable in such a case.  \
As with the "inputs" list, the variable names should be strings and the variable types should be ADVarType objects.

operation: a function implementing the intended operation within python.  This function must take, as arguments, the \
inputs defined by the "inputs" list, and must return a tuple of numpy.ndarray objects representing the outputs defined \
by the "outputs" list.  The argument names must match the names given in the "inputs" list and should appear in the \
same order as defined in the "inputs" list.  The outputs must be returned with the same order defined by the \
"outputs" list.
IMPORTANT: all outputs must be returned in a tuple of non-scalar numpy.ndarray objects, even scalar outputs.  Any \
scalar outputs should be contained in a 1-dimensional, length 1 numpy.ndarray.  For example, an operation that has a \
single scalar output should return "(numpy.array([scalar_output]), )".

variants: (optional) a list of sets of input variable names defining the variations of the function to test.  Each \
set of input variable names defines a test case where the given inputs are assumed to be differentiable, all other \
inputs are assumed to be constants.  If left undefined, then the power set of variable names from "inputs" is used.
"""

    # Uses "/*UNQ_T1F_<method acronym>_##*/" as the unique identifier comment for debugging purposes.

    __hc__ = 1e-20
    __h__ = 1e-5  # these values give the most accurate performance

    __test_abs_err__ = {}

    operation_name: str
    operation: operation_type
    inputs: io_ad_type  # [(<var_name>, <var_type>), ...]
    outputs: io_ad_type  # [(<var_name>, <var_type>), ...]
    variants: variants_type

    def __init__(self,
                 operation_name: str,
                 inputs: io_ad_type,
                 outputs: io_ad_type,
                 operation: operation_type,
                 variants: variants_input_type = None,
                 ):
        self.operation_name = operation_name
        self.operation = operation
        self.inputs = inputs
        self.outputs = outputs

        # Make sure inputs and outputs all have unique names
        var_names = {var_name for var_name, var_type in itertools.chain(self.inputs, self.outputs)}
        if len(var_names) != (len(self.inputs) + len(self.outputs)):
            raise VariableNameOverloadError(f'Overlap detected in variable names for\n'
                                            f'{self.__class__.__name__}({self.operation_name}, ...):\n'
                                            f'\tinputs = {self.inputs}\n'
                                            f'\toutputs = {self.outputs}\n')

        # Check and/or assign functional variants
        if variants is not None:
            # noinspection PyTypeChecker
            self.variants = sorted(map(frozenset, variants))
            expected_variants = self.__expected_variants__()
            frozen_variants = set(self.variants)
            if bool(expected_variants.symmetric_difference(frozen_variants)):
                w = UserWarning(f"\nWARNING: Incomplete or unexpected set of function variants prescribed to:\n"
                                f'\t\t{self.__class__.__name__}("{self.operation_name}", ...)\n'
                                f"\tGiven: {tuple(map(set, frozen_variants))}\n"
                                f"\tExpected: {tuple(map(set, expected_variants))}\n")
                warnings.warn(w)
                pass
            pass
        else:
            # "variants" now falls back on __expected_variants__
            self.variants = sorted(self.__expected_variants__())
            pass

        # Inspect operation to check for agreement between inputs and outputs
        self.__inspect_operation__()

        pass

    def __inspect_operation__(self) -> None:
        """Make sure the operation signature matches that of the inputs and outputs."""
        operation_arg_spec = inspect.getfullargspec(self.operation)

        # Deal with variadic arguments
        if operation_arg_spec.varargs is not None:
            # TODO
            w = FutureWarning(f'from TestADFunction({self.operation_name}, ...)'
                              f'\n\tTestADFunction.__inspect_operation__ is not yet capable of inspecting variadic'
                              f'\n\targuments in the "operation" function.  Dependable inspection and agreement'
                              f'\n\tchecking is unlikely.')
            warnings.warn(w)
            pass

        # Deal with variadic keyword arguments
        if operation_arg_spec.varkw is not None:
            # TODO
            w = FutureWarning(f'from TestADFunction({self.operation_name}, ...)'
                              f'\n\tTestADFunction.__inspect_operation__ is not yet capable of inspecting variadic'
                              f'\n\tkeyword arguments in the "operation" function.  Dependable inspection and'
                              f'\n\tagreement checking is unlikely.')
            warnings.warn(w)
            pass

        # Check to see if the function is annotated
        if operation_arg_spec.annotations:
            # function is annotated
            # TODO
            pass
        else:
            # function is not annotated
            # TODO
            pass
        pass

    def __expected_variants__(self) -> Set[FrozenSet[str]]:
        input_names = [var_name for var_name, var_type in self.inputs]
        input_name_power_set = set(
            map(frozenset, itertools.chain.from_iterable(itertools.combinations(input_names, r)
                                                         for r in range(len(input_names) + 1)))
        )
        return input_name_power_set

    @property
    def output_names(self) -> Tuple[str, ...]:
        return next(zip(*self.outputs))

    def evaluate_at(self, input_data: input_data_type):
        kwargs = {var_name: input_data[var_name][0]
                  for var_name, var_type in self.inputs}
        return self.operation(**kwargs)

    def evaluate_output_derivatives_at(self, input_data: input_data_type,
                                       non_constant_inputs: non_constant_inputs_type):
        kwargs = {
            var_name: (
                input_data[var_name][0] + (1j * self.__hc__ * input_data[f"{var_name}b"][0])
                if var_name in non_constant_inputs else
                input_data[var_name][0]
            )
            for var_name, var_type in self.inputs
        }
        return [numpy.imag(x) / self.__hc__ for x in self.operation(**kwargs)]

    def evaluate_input_derivatives_at(self, input_data: input_data_type,
                                      non_constant_inputs: non_constant_inputs_type):
        output_names = self.output_names
        dimension_vectors = {
            var_name: [numpy.reshape(dim, numpy.shape(input_data[var_name][0]))
                       for dim in numpy.identity(numpy.size(input_data[var_name][0]))]
            for var_name, var_type in self.inputs
            if var_name in non_constant_inputs
        }
        partial_derivatives = {
            candidate_var_name: [
                {output_name: numpy.imag(x) / self.__hc__
                 for output_name, x in
                 zip(output_names, self.operation(**{
                     var_name: (
                         input_data[var_name][0] + (1j * self.__hc__ * dimension_vector)
                         if var_name == candidate_var_name else
                         input_data[var_name][0]
                     )
                     for var_name, var_type in self.inputs
                 }))}
                for dimension_vector in dimension_vectors[candidate_var_name]
            ]
            for candidate_var_name, candidate_type in self.inputs
            if candidate_var_name in non_constant_inputs
        }

        results = {
            candidate_var_name: numpy.reshape([
                sum(numpy.sum(numpy.multiply(input_data[f"{output_name}b"][0], derivative))
                    for output_name, derivative in partial.items())
                for partial in partials
            ], numpy.shape(input_data[candidate_var_name][0]))
            for candidate_var_name, partials in partial_derivatives.items()
        }

        return [results[var_name]
                for var_name, var_type in self.inputs
                if var_name in non_constant_inputs]

    @property
    def datatypes(self) -> Set[ADVarType]:
        # datatypes = set(list(zip(*self.inputs))[1]).union(list(zip(*self.outputs))[1])
        datatypes = self.input_types
        datatypes.update(self.output_types)
        return datatypes

    @property
    def input_types(self) -> Set[ADVarType]:
        return {var_type for var_name, var_type in self.inputs}

    @property
    def output_types(self) -> Set[ADVarType]:
        return {var_type for var_name, var_type in self.outputs}

    def _test_variant_name_(self, non_constant_inputs: non_constant_inputs_type) -> str:
        modifier = ''.join(var_type.shortname if var_name in non_constant_inputs else var_type.parent.shortname
                           for var_name, var_type in self.inputs
                           # if var_name in non_constant_inputs
                           )
        return f"{self.operation_name}_{modifier}"

    @staticmethod
    def _all_declaration_format_(string_format: str, io_vars: io_ad_type, var_type: ADVarType, **kwargs) -> List[str]:
        defaults = dict(zip(['default_value', 'default_bvalue'], var_type.initialization_defaults))
        return [string_format.format(var_name=var_name, **defaults, **kwargs)
                for var_name, io_var_type in io_vars
                if io_var_type == var_type]  # This might need modification

    @staticmethod
    def _non_constant_declaration_format_(string_format: str, io_vars: io_ad_type, var_type: ADVarType,
                                          non_constant_inputs: non_constant_inputs_type, **kwargs) -> List[str]:
        defaults = dict(zip(['default_value', 'default_bvalue'], var_type.initialization_defaults))
        return [string_format.format(var_name=var_name, **defaults, **kwargs)
                for var_name, io_var_type in io_vars
                if (io_var_type == var_type and  # This might need modification
                    var_name in non_constant_inputs)]

    def _evaluation_signature_(self, non_constant_inputs: non_constant_inputs_type) -> str:
        return ', '.join(
            [  # Inputs
                f"{var_name}_ad" if (var_name in non_constant_inputs) else f"{var_name}"
                for var_name, var_type in self.inputs
            ] +
            [  # Outputs
                f"{var_name}_ad" if non_constant_inputs else f"{var_name}"
                # ^ If any of the inputs are non-constant, then all the outputs are non-constant
                for var_name, var_type in self.outputs
            ]
        )

    @staticmethod
    def _test_f_(test_variant_name: str, test_type: str,
                 declarations_list: List[str], evaluations_list: List[str], comparisons_list: List[str]) -> str:
        separator_string = '\n    '
        declarations = separator_string.join(declarations_list)
        evaluations = separator_string.join(evaluations_list)
        comparisons = separator_string.join(comparisons_list)
        return f"""
TEST_F({test_variant_name}, {test_type}) {{
    // Declarations and Initializations:
    {declarations}
    // Evaluations:
    {evaluations}
    // Comparisons:
    {comparisons}
}}
"""

    def _test_f_passive_(self, test_variant_name: str, non_constant_inputs: non_constant_inputs_type) -> str:
        non_constant_input_types = {var_type for var_name, var_type in self.inputs if var_name in non_constant_inputs}
        declarations_list = (
                [  # initialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard inputs
                        self._all_declaration_format_(
                            "{var_name}({var_name}_data)", self.inputs, var_type
                        )
                    ) +
                    ";  /*UNQ_T1F_TFP_01*/"
                    for var_type in self.input_types
                ] +
                [  # uninitialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard outputs
                        self._all_declaration_format_(
                            "{var_name}", self.outputs, var_type
                        ) +
                        # "b" outputs (FYI: if any of the inputs are non-constant, then all outputs are non-constant)
                        (self._all_declaration_format_(
                            "{var_name}b", self.outputs, var_type
                        ) * bool(non_constant_inputs)) +
                        # "b" inputs (for the non-constant inputs)
                        self._non_constant_declaration_format_(
                            "{var_name}b", self.inputs, var_type, non_constant_inputs
                        )
                    ) +
                    ";  /*UNQ_T1F_TFP_02*/"
                    for var_type in self.datatypes
                    # Make sure this declaration isn't empty
                    if ((var_type in self.output_types or
                         var_type in non_constant_input_types) and
                        (var_type.requires_parent_declaration or
                         not bool(non_constant_inputs)))
                ] +
                [  # dependent declarations (A.K.A. AD initializations)
                    f"{var_type.alias} " +
                    ', '.join(
                        self._non_constant_declaration_format_(
                            "{var_name}_ad({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_ad({var_name}, {default_bvalue})",
                            self.inputs, var_type, non_constant_inputs
                        ) +
                        (self._all_declaration_format_(
                            "{var_name}_ad({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_ad({default_value}, {default_bvalue})",
                            self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T1F_TFP_03*/"
                    for var_type in self.datatypes
                    # Make sure this declaration isn't empty
                    if ((var_type in self.output_types and bool(non_constant_inputs)) or
                        var_type in non_constant_input_types)
                ]
        )

        evaluations_list = [
            # f'{auto_typing}A2D::{self.operation_name}({self._evaluation_signature_(non_constant_inputs)})'
            f'A2D::{self.operation_name}({self._evaluation_signature_(non_constant_inputs)})'
            f';  /*UNQ_T1F_TFP_04*/'
        ]

        comparisons_list = (  # This depends on type of test (e.g. passive, forward, etc...)
            [
                f"{var_type.test}({var_name}_ad{var_type.get_value}, {var_name}_out);  /*UNQ_T1F_TFP_05*/"
                for var_name, var_type in self.outputs
            ] if bool(non_constant_inputs) else
            [
                f"{var_type.test}({var_name}, {var_name}_out);  /*UNQ_T1F_TFP_06*/"
                for var_name, var_type in self.outputs
            ]
        )
        return self._test_f_(test_variant_name, "passive", declarations_list, evaluations_list, comparisons_list)

    def _test_f_forward_(self, test_variant_name: str, non_constant_inputs: non_constant_inputs_type) -> str:
        assert bool(non_constant_inputs)  # this might need to change later
        non_constant_input_types = {var_type for var_name, var_type in self.inputs if var_name in non_constant_inputs}
        declarations_list = (
                [  # initialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard inputs
                        self._all_declaration_format_(
                            "{var_name}({var_name}_data)", self.inputs, var_type
                        ) +
                        # "b" inputs (for the non-constant inputs)
                        self._non_constant_declaration_format_(
                            "{var_name}b({var_name}b_data)", self.inputs, var_type, non_constant_inputs
                        )
                    ) +
                    ";  /*UNQ_T1F_TFF_01*/"
                    for var_type in self.input_types
                ] +
                [  # uninitialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard outputs
                        self._all_declaration_format_(
                            "{var_name}", self.outputs, var_type
                        ) +
                        # "b" outputs (FYI: if any of the inputs are non-constant, then all outputs are non-constant)
                        (self._all_declaration_format_(
                            "{var_name}b", self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T1F_TFF_02*/"
                    for var_type in self.output_types
                    if var_type.requires_parent_declaration
                ] +
                [  # dependent declarations (A.K.A. AD initializations)
                    f"{var_type.alias} " +
                    ', '.join(
                        self._non_constant_declaration_format_(
                            "{var_name}_ad({var_name}, {var_name}b)",
                            self.inputs, var_type, non_constant_inputs
                        ) +
                        (self._all_declaration_format_(
                            "{var_name}_ad({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_ad({default_value}, {default_bvalue})",
                            self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T1F_TFF_03*/"
                    for var_type in self.datatypes
                    # Make sure this declaration isn't empty
                    if var_type in self.output_types or var_type in non_constant_input_types
                ]
        )

        evaluations_list = [
            f'auto expr = A2D::{self.operation_name}({self._evaluation_signature_(non_constant_inputs)})'
            f';  /*UNQ_T1F_TFF_04*/',
            f'expr.forward();  /*UNQ_T1F_TFF_05*/',
        ]

        comparisons_list = (  # This depends on the type of test (e.g. passive, forward, etc...)
                [
                    f"{var_type.test}({var_name}_ad{var_type.get_value}, {var_name}_out);  /*UNQ_T1F_TFF_06*/"
                    for var_name, var_type in self.outputs
                ] +
                [
                    f"{var_type.test}({var_name}_ad{var_type.get_bvalue}, {var_name}b_out);  /*UNQ_T1F_TFF_07*/"
                    for var_name, var_type in self.outputs
                ]
        )
        return self._test_f_(test_variant_name, "forward", declarations_list, evaluations_list, comparisons_list)

    def _test_f_reverse_(self, test_variant_name: str, non_constant_inputs: non_constant_inputs_type) -> str:
        assert bool(non_constant_inputs)  # this might need to change later
        non_constant_input_types = {var_type for var_name, var_type in self.inputs if var_name in non_constant_inputs}
        declarations_list = (
                [  # initialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard inputs
                        self._all_declaration_format_(
                            "{var_name}({var_name}_data)", self.inputs, var_type
                        ) +
                        # "b" outputs (FYI: if any of the inputs are non-constant, then all outputs are non-constant)
                        (self._all_declaration_format_(
                            "{var_name}b({var_name}b_data)", self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T1F_TFR_01*/"
                    for var_type in self.datatypes
                ] +
                [  # uninitialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard outputs
                        self._all_declaration_format_(
                            "{var_name}", self.outputs, var_type
                        ) +
                        # "b" inputs (for the non-constant inputs)
                        self._non_constant_declaration_format_(
                            "{var_name}b", self.inputs, var_type, non_constant_inputs
                        )
                    ) +
                    ";  /*UNQ_T1F_TFR_02*/"
                    for var_type in self.datatypes
                    if (var_type.requires_parent_declaration and
                        (var_type in non_constant_input_types or
                         var_type in self.output_types))
                ] +
                [  # dependent declarations (A.K.A. AD initializations)
                    f"{var_type.alias} " +
                    ', '.join(
                        self._non_constant_declaration_format_(
                            "{var_name}_ad({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_ad({var_name}, {default_bvalue})",
                            self.inputs, var_type, non_constant_inputs
                        ) +
                        (self._all_declaration_format_(
                            "{var_name}_ad({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_ad({default_value}, {var_name}b)",
                            self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T1F_TFR_03*/"
                    for var_type in self.datatypes
                    # Make sure this declaration isn't empty
                    if var_type in self.output_types or var_type in non_constant_input_types
                ]
        )

        evaluations_list = [
            f'auto expr = A2D::{self.operation_name}({self._evaluation_signature_(non_constant_inputs)})'
            f';  /*UNQ_T1F_TFR_04*/',
            f'expr.reverse();  /*UNQ_T1F_TFR_05*/',
        ]

        comparisons_list = (  # This depends on type of test (e.g. passive, forward, etc...)
                [
                    f"{var_type.test}({var_name}_ad{var_type.get_value}, {var_name}_out);  /*UNQ_T1F_TFR_06*/"
                    for var_name, var_type in self.outputs
                ] +
                [
                    f"{var_type.test}({var_name}_ad{var_type.get_bvalue}, {var_name}b_out);  /*UNQ_T1F_TFR_07*/"
                    for var_name, var_type in self.inputs
                    if var_name in non_constant_inputs
                ]
        )
        return self._test_f_(test_variant_name, "reverse", declarations_list, evaluations_list, comparisons_list)

    def full_test_variant(self, input_data: input_data_type,
                          non_constant_inputs: non_constant_inputs_type) -> str:
        test_variant_name = self._test_variant_name_(non_constant_inputs)
        separator_string = '\n    '

        # generate outputs
        if non_constant_inputs:
            nc_inputs = [(var_name, var_type)
                         for var_name, var_type in self.inputs
                         if var_name in non_constant_inputs]
            # Include derivatives
            results = (
                    [
                        # calculated nominal values
                        var_type.initialize_assignment_format.format(var_name=f"{var_name}_out", *value)
                        for (var_name, var_type), value in zip(self.outputs, self.evaluate_at(input_data))
                    ] +
                    [
                        # calculated output derivatives
                        var_type.initialize_assignment_format.format(var_name=f"{var_name}b_out", *value)
                        for (var_name, var_type), value in
                        zip(self.outputs, self.evaluate_output_derivatives_at(input_data, non_constant_inputs))
                    ] +
                    [
                        # calculated input derivatives
                        var_type.initialize_assignment_format.format(var_name=f"{var_name}b_out", *value)
                        for (var_name, var_type), value in
                        zip(nc_inputs, self.evaluate_input_derivatives_at(input_data, non_constant_inputs))
                    ]
            )
        else:
            # No derivatives necessary
            results = (
                [
                    # calculated nominal values
                    var_type.initialize_assignment_format.format(var_name=f"{var_name}_out", *value)
                    for (var_name, var_type), value in zip(self.outputs, self.evaluate_at(input_data))
                ]
            )

        test_string = f"""
class {test_variant_name} : public {self.operation_name} {{
protected:
    // Results
    {separator_string.join(results)}
}};  /*UNQ_T1F_FTV_01*/
{self._test_f_passive_(test_variant_name, non_constant_inputs)}"""

        if non_constant_inputs:
            return test_string + \
                   self._test_f_forward_(test_variant_name, non_constant_inputs) + \
                   self._test_f_reverse_(test_variant_name, non_constant_inputs)
        return test_string

    def tests(self, input_data: input_data_type) -> str:
        tests = (self.full_test_variant(input_data, non_constant_inputs)
                 for non_constant_inputs in self.variants)
        return ''.join(tests)

    pass


class TestA2DFunction(TestADFunction):
    """\
TestA2DFunction: A representation of an overloaded operation to be tested.

operation_name: the name of the overloaded operation within the A2D C++ code as a string using double quotations.  \
The operation should return an instance of the resulting A2D class, specifically an "A2D" class; all outputs and \
modified objects must be passed by reference to the operation.

inputs: a list of variable-name, variable-type pairs that define the inputs to the overloaded operation in the MOST \
differentiable state, that is, when no input is considered constant.  The variable names should be strings, and the \
variable types should be A2DVarType objects.

outputs: a list of variable-name, variable-type pairs that define the outputs to the overloaded operation in the MOST \
differentiable state.  For most if not all circumstances, all outputs should be differentiable in such a case.  \
As with the "inputs" list, the variable names should be strings and the variable types should be A2DVarType objects.

operation: a function implementing the intended operation within python.  This function must take, as arguments, the \
inputs defined by the "inputs" list, and must return a tuple of numpy.ndarray objects representing the outputs defined \
by the "outputs" list.  The argument names must match the names given in the "inputs" list and should appear in the \
same order as defined in the "inputs" list.  The outputs must be returned with the same order defined by the \
"outputs" list.
IMPORTANT: all outputs must be returned in a tuple of non-scalar numpy.ndarray objects, even scalar outputs.  Any \
scalar outputs should be contained in a 1-dimensional, length 1 numpy.ndarray.  For example, an operation that has a \
single scalar output should return "(numpy.array([scalar_output]), )".

variants: (optional) a list of sets of input variable names defining the variations of the function to test.  Each \
set of input variable names defines a test case where the given inputs are assumed to be differentiable, all other \
inputs are assumed to be constants.  If left undefined, then the power set of variable names from "inputs" is used.
"""

    # Uses "/*UNQ_T2F_<method acronym>_##*/" as the unique identifier comment for debugging purposes.

    operation_name: str
    operation: operation_type
    inputs: io_a2d_type  # [(<var_name>, <var_type>), ...]
    outputs: io_a2d_type  # [(<var_name>, <var_type>), ...]
    variants: variants_type
    N: int

    # noinspection PyPep8Naming
    def __init__(self,
                 operation_name: str,
                 inputs: io_a2d_type,
                 outputs: io_a2d_type,
                 operation: operation_type,
                 variants: variants_input_type = None,
                 ):
        TestADFunction.__init__(self, operation_name, inputs, outputs, operation, variants)
        datatypes = self.datatypes.copy()
        N = datatypes.pop().N
        if any(var_type.N != N for var_type in datatypes):
            N_values = {N}.union(var_type.N for var_type in datatypes)
            raise ValueError(f'from TestA2DFunction({operation_name}, ...)\n'
                             f'All values of "N" must be equal for all input and output A2DVarTypes.\n'
                             f'Instead, there were {len(N_values)} different values for "N": {N_values}')
        self.N = N
        pass

    @property
    def datatypes(self) -> Set[A2DVarType]:
        datatypes = set(list(zip(*self.inputs))[1]).union(list(zip(*self.outputs))[1])
        return datatypes

    @property
    def input_types(self) -> Set[A2DVarType]:
        return {var_type for var_name, var_type in self.inputs}

    @property
    def output_types(self) -> Set[A2DVarType]:
        return {var_type for var_name, var_type in self.outputs}

    def __expected_variants__(self):
        input_name_power_set = TestADFunction.__expected_variants__(self)
        return input_name_power_set.difference({frozenset()})

    def evaluate_output_derivatives_at(self, input_data: input_data_type,
                                       non_constant_inputs: non_constant_inputs_type):
        results = []
        for i in range(self.N):
            kwargs = {
                var_name: (
                    input_data[var_name][0] + (1j * self.__hc__ * input_data[f"{var_name}p{i}"][0])
                    if var_name in non_constant_inputs else
                    input_data[var_name][0]
                )
                for var_name, var_type in self.inputs
            }
            results.append([numpy.imag(x) / self.__hc__ for x in self.operation(**kwargs)])
        return results

    def __flat_operation__(self, input_vector: numpy.ndarray,
                           input_shapes: List[Tuple[str, Tuple[int, ...]]],
                           input_sizes: List[int]):
        split_indices = [0] + list(itertools.accumulate(input_sizes))
        kwargs = {
            var_name: numpy.reshape(input_vector[s1:s2], var_shape)
            for (var_name, var_shape), (s1, s2) in
            zip(input_shapes, zip(split_indices, split_indices[1:]))
        }
        return numpy.concatenate([res.flatten() for res in self.operation(**kwargs)])

    # noinspection PyPep8Naming
    def evaluate_hessian_derivatives_at(self, input_data: input_data_type,
                                        non_constant_inputs: non_constant_inputs_type
                                        ) -> List[List[Tuple[str, A2DVarType, numpy.ndarray]]]:
        # determine input shapes
        input_shapes = [(var_name, numpy.shape(input_data[var_name][0]))
                        for var_name, var_type in self.inputs]
        input_sizes = [numpy.size(input_data[var_name][0])
                       for var_name, var_type in self.inputs]
        input_size = sum(input_sizes)
        # determine non constant input shapes
        nc_inputs = [(var_name, var_type)
                     for var_name, var_type in self.inputs
                     if var_name in non_constant_inputs]
        nc_input_sizes = [var_size
                          for var_size, (var_name, var_type) in
                          zip(input_sizes, self.inputs)
                          if var_name in non_constant_inputs]

        # nominal inputs and delta matrix
        nominal_input_vector = numpy.concatenate([input_data[var_name][0].flatten()
                                                  for var_name, var_type in self.inputs])
        non_constant_indices = numpy.concatenate([[var_name in non_constant_inputs] * size
                                                  for (var_name, var_type), size in
                                                  zip(self.inputs, input_sizes)])
        delta_mat = numpy.identity(input_size)[non_constant_indices]

        # load and flatten b values
        flat_Cb = numpy.concatenate([input_data[f"{var_name}b"][0].flatten()
                                     for var_name, var_type in self.outputs])
        # compute jacobian matrix
        flat_dCdV = numpy.array([
            numpy.imag(self.__flat_operation__(nominal_input_vector + delta, input_shapes, input_sizes)
                       ) / self.__hc__
            for delta in self.__hc__ * 1j * delta_mat
        ])
        # compute hessian tensor
        flat_d2CdV2 = numpy.array([
            [
                numpy.imag(
                    self.__flat_operation__(nominal_input_vector + delta_a + delta_b, input_shapes, input_sizes) -
                    self.__flat_operation__(nominal_input_vector - delta_a + delta_b, input_shapes, input_sizes)
                ) / (2 * self.__hc__ * self.__h__)
                for delta_b in self.__hc__ * 1j * delta_mat
            ]
            for delta_a in self.__h__ * delta_mat
        ])
        results = []
        for i in range(self.N):
            flat_Vp = numpy.concatenate([input_data[f"{var_name}p{i}"][0].flatten()
                                         for var_name, var_type in nc_inputs])
            flat_Ch = numpy.concatenate([input_data[f"{var_name}h{i}"][0].flatten()
                                         for var_name, var_type in self.outputs])
            flat_Vh = numpy.dot(numpy.dot(flat_d2CdV2, flat_Cb), flat_Vp) + numpy.dot(flat_dCdV, flat_Ch)

            # separate out the results
            split_indices = [0] + list(itertools.accumulate(nc_input_sizes))
            Vh = [flat_Vh[s1:s2]
                  for s1, s2 in
                  zip(split_indices, split_indices[1:])]

            results.append([(var_name, var_type, value)
                            for (var_name, var_type), value in
                            zip(nc_inputs, Vh)])
        return results

    @staticmethod
    def _construct_derivative_loop_(derivative_loops_list: List[List[str]]) -> str:
        indent_string = '    '
        max_rank = len(derivative_loops_list) - 1
        index_list = []
        derivative_loop_string = "/*UNQ_T2F_CDL_01*/"
        for rank, derivative_assignments in enumerate(derivative_loops_list[:-1]):
            separator_string = '\n' + (indent_string * rank)
            if derivative_assignments:
                derivative_loop_string += separator_string
                derivative_loop_string += separator_string.join(
                    derivative_assignment.format(**{f"rank{rank}": ('(' + ', '.join(index_list) + ')') * bool(rank)})
                    for derivative_assignment in derivative_assignments
                )
            derivative_loop_string += (separator_string +
                                       f'for (int ii_{rank} = 0; ii_{rank} < 3; ii_{rank}++) {{  /*UNQ_T2F_CDL_02*/')
            index_list.append(f'ii_{rank}')
        separator_string = '\n' + (indent_string * max_rank)
        derivative_loop_string += separator_string
        derivative_loop_string += separator_string.join(
            derivative_assignment.format(**{f"rank{max_rank}": ('(' + ', '.join(index_list) + ')') * bool(max_rank)})
            for derivative_assignment in derivative_loops_list[-1]
        )

        # close brackets:
        for rank in reversed(range(max_rank)):
            derivative_loop_string += '\n' + (indent_string * rank) + '}'
        return derivative_loop_string

    def _evaluation_signature_(self, non_constant_inputs: non_constant_inputs_type):
        return ', '.join(
            [  # Inputs
                f"{var_name}_a2d" if (var_name in non_constant_inputs) else f"{var_name}"
                for var_name, var_type in self.inputs
            ] +
            [  # Outputs
                f"{var_name}_a2d" if non_constant_inputs else f"{var_name}"
                # ^ If any of the inputs are non-constant, then all the outputs are non-constant
                for var_name, var_type in self.outputs
            ]
        )

    # noinspection PyMethodOverriding
    @staticmethod
    def _test_f_(test_variant_name: str, test_type: str,
                 declarations_list: List[str], derivative_loop: str,
                 evaluations_list: List[str], comparisons_list: List[str]) -> str:
        separator_string = '\n    '
        declarations = separator_string.join(declarations_list)
        # derivative_loops = separator_string.join(
        #     [derivative_loop.replace('\n', separator_string) for derivative_loop in derivative_loops_list]
        # )
        derivative_loop = derivative_loop.replace('\n', separator_string)
        evaluations = separator_string.join(evaluations_list)
        comparisons = separator_string.join(comparisons_list)
        return f"""
TEST_F({test_variant_name}, {test_type}) {{
    // Declarations and Initializations:
    {declarations}
    // Set Derivative Values:
    {derivative_loop}
    // Evaluations:
    {evaluations}
    // Comparisons:
    {comparisons}
}}
"""

    def _test_f_passive_(self, test_variant_name: str, non_constant_inputs: non_constant_inputs_type) -> str:
        non_constant_input_types = {var_type for var_name, var_type in self.inputs if var_name in non_constant_inputs}
        declarations_list = (
                [  # initialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard inputs
                        self._all_declaration_format_(
                            "{var_name}({var_name}_data)", self.inputs, var_type
                        )
                    ) +
                    ";  /*UNQ_T2F_TFP_01*/"
                    for var_type in self.input_types
                ] +
                [  # uninitialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard outputs
                        self._all_declaration_format_(
                            "{var_name}", self.outputs, var_type
                        ) +
                        # "b" outputs (FYI: if any of the inputs are non-constant, then all outputs are non-constant)
                        (self._all_declaration_format_(
                            "{var_name}b", self.outputs, var_type
                        ) * bool(non_constant_inputs)) +
                        # "b" inputs (for the non-constant inputs)
                        self._non_constant_declaration_format_(
                            "{var_name}b", self.inputs, var_type, non_constant_inputs
                        )
                    ) +
                    ";  /*UNQ_T2F_TFP_02*/"
                    for var_type in self.datatypes
                    # Make sure this declaration isn't empty
                    if ((var_type in self.output_types or
                         var_type in non_constant_input_types) and
                        (var_type.requires_parent_declaration or
                         not bool(non_constant_inputs)))
                ] +
                [  # dependent declarations (A.K.A. A2D initializations)
                    f"{var_type.alias} " +
                    ', '.join(
                        self._non_constant_declaration_format_(
                            "{var_name}_a2d({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_a2d({var_name}, {default_bvalue})",
                            self.inputs, var_type, non_constant_inputs
                        ) +
                        (self._all_declaration_format_(
                            "{var_name}_a2d({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_a2d({default_value}, {default_bvalue})",
                            self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T2F_TFP_03*/"
                    for var_type in self.datatypes
                    # Make sure this declaration isn't empty
                    if ((var_type in self.output_types and bool(non_constant_inputs)) or
                        var_type in non_constant_input_types)
                ]
        )

        derivative_loop = '    /*None for "passive" tests*/'

        evaluations_list = [
            f'A2D::{self.operation_name}({self._evaluation_signature_(non_constant_inputs)})'
            f';  /*UNQ_T2F_TFP_04*/'
        ]

        comparisons_list = (  # This depends on type of test (e.g. passive, forward, etc...)
            [
                f"{var_type.test}({var_name}_a2d{var_type.get_value}, {var_name}_out);  /*UNQ_T2F_TFP_05*/"
                for var_name, var_type in self.outputs
            ] if bool(non_constant_inputs) else
            [
                f"{var_type.test}({var_name}, {var_name}_out);  /*UNQ_T2F_TFP_06*/"
                for var_name, var_type in self.outputs
            ]
        )
        return self._test_f_(test_variant_name, "passive", declarations_list,
                             derivative_loop, evaluations_list, comparisons_list)

    def _test_f_reverse_(self, test_variant_name: str, non_constant_inputs: non_constant_inputs_type) -> str:
        assert bool(non_constant_inputs)  # this might need to change later
        non_constant_input_types = {var_type for var_name, var_type in self.inputs if var_name in non_constant_inputs}
        declarations_list = (
                [  # initialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard inputs
                        self._all_declaration_format_(
                            "{var_name}({var_name}_data)", self.inputs, var_type
                        ) +
                        # "b" outputs (FYI: if any of the inputs are non-constant, then all outputs are non-constant)
                        (self._all_declaration_format_(
                            "{var_name}b({var_name}b_data)", self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T2F_TFR_01*/"
                    for var_type in self.datatypes
                ] +
                [  # uninitialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(
                        # standard outputs
                        self._all_declaration_format_(
                            "{var_name}", self.outputs, var_type
                        ) +
                        # "b" inputs (for the non-constant inputs)
                        self._non_constant_declaration_format_(
                            "{var_name}b", self.inputs, var_type, non_constant_inputs
                        )
                    ) +
                    ";  /*UNQ_T2F_TFR_02*/"
                    for var_type in self.datatypes
                    if (var_type.requires_parent_declaration and
                        (var_type in non_constant_input_types or
                         var_type in self.output_types))
                ] +
                [  # dependent declarations (A.K.A. A2D initializations)
                    f"{var_type.alias} " +
                    ', '.join(
                        self._non_constant_declaration_format_(
                            "{var_name}_a2d({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_a2d({var_name}, {default_bvalue})",
                            self.inputs, var_type, non_constant_inputs
                        ) +
                        (self._all_declaration_format_(
                            "{var_name}_a2d({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_a2d({default_value}, {var_name}b)",
                            self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T2F_TFR_03*/"
                    for var_type in self.datatypes
                    # Make sure this declaration isn't empty
                    if var_type in self.output_types or var_type in non_constant_input_types
                ]
        )

        derivative_loop = '    /*None for "reverse" tests*/'

        evaluations_list = [
            f'auto expr = A2D::{self.operation_name}({self._evaluation_signature_(non_constant_inputs)})'
            f';  /*UNQ_T2F_TFR_04*/',
            f'expr.reverse();  /*UNQ_T2F_TFR_05*/',
        ]

        comparisons_list = (  # This depends on type of test (e.g. passive, forward, etc...)
                [
                    f"{var_type.test}({var_name}_a2d{var_type.get_value}, {var_name}_out);  /*UNQ_T2F_TFR_06*/"
                    for var_name, var_type in self.outputs
                ] +
                [
                    f"{var_type.test}({var_name}_a2d{var_type.get_bvalue}, {var_name}b_out);  /*UNQ_T2F_TFR_07*/"
                    for var_name, var_type in self.inputs
                    if var_name in non_constant_inputs
                ]
        )
        return self._test_f_(test_variant_name, "reverse", declarations_list,
                             derivative_loop, evaluations_list, comparisons_list)

    def _test_f_h_forward_(self, test_variant_name: str, non_constant_inputs: non_constant_inputs_type) -> str:
        assert bool(non_constant_inputs)  # this might need to change later
        non_constant_input_types = {var_type for var_name, var_type in self.inputs if var_name in non_constant_inputs}
        declarations_list = (
                [  # initialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(itertools.chain(
                        # standard inputs
                        self._all_declaration_format_(
                            "{var_name}({var_name}_data)", self.inputs, var_type
                        ),
                        # "p" inputs (for the non-constant inputs)
                        *(self._non_constant_declaration_format_(
                            "{var_name}p{i}({var_name}p{i}_data)", self.inputs, var_type, non_constant_inputs, i=i
                        ) for i in range(self.N)),
                        # "b" outputs (FYI: if any of the inputs are non-constant, then all outputs are non-constant)
                        (self._all_declaration_format_(
                            "{var_name}b({var_name}b_data)", self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    )) +
                    ";  /*UNQ_T2F_TFHF_01*/"
                    for var_type in self.datatypes
                ] +
                [  # uninitialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(itertools.chain(
                        # standard outputs
                        self._all_declaration_format_(
                            "{var_name}", self.outputs, var_type
                        ),
                        # "b" inputs (for the non-constant inputs)
                        self._non_constant_declaration_format_(
                            "{var_name}b", self.inputs, var_type, non_constant_inputs
                        )
                    )) +
                    ";  /*UNQ_T2F_TFHF_02*/"
                    # for var_type in self.output_types
                    # if var_type.requires_parent_declaration
                    for var_type in self.datatypes
                    if (var_type.requires_parent_declaration and
                        (var_type in non_constant_input_types or
                         var_type in self.output_types))
                ] +
                [  # dependent declarations (A.K.A. A2D initializations)
                    f"{var_type.alias} " +
                    ', '.join(
                        self._non_constant_declaration_format_(
                            "{var_name}_a2d({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_a2d({var_name}, {default_bvalue})",
                            self.inputs, var_type, non_constant_inputs
                        ) +
                        (self._all_declaration_format_(
                            "{var_name}_a2d({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_a2d({default_value}, {var_name}b)",
                            self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T2F_TFHF_03*/"
                    for var_type in self.datatypes
                    # Make sure this declaration isn't empty
                    if var_type in self.output_types or var_type in non_constant_input_types
                ]
        )

        derivative_loops_list = [
            [  # p value assignments
                f"{var_name}_a2d{var_type.get_pvalue.format(i)}{{rank{rank}}} = {var_name}p{i}{{rank{rank}}}"
                f";  /*UNQ_T2F_TFHF_05*/"
                for i in range(self.N)
                for var_name, var_type in self.inputs
                if (var_name in non_constant_inputs and
                    var_type.rank == rank)
            ] +
            [  # h value assignments
            ]
            for rank in range(max(inp.rank for inp in non_constant_input_types) + 1)
        ]
        derivative_loop = self._construct_derivative_loop_(derivative_loops_list)

        evaluations_list = [
            f'auto expr = A2D::{self.operation_name}({self._evaluation_signature_(non_constant_inputs)})'
            f';  /*UNQ_T2F_TFHF_06*/',
            f'expr.reverse();  /*UNQ_T2F_TFHF_07*/',
            f'expr.hforward();  /*UNQ_T2F_TFHF_08*/',
        ]

        comparisons_list = list(itertools.chain(  # This depends on the type of test (e.g. passive, forward, etc...)
            [  # value check
                f"{var_type.test}({var_name}_a2d{var_type.get_value}, {var_name}_out);  /*UNQ_T2F_TFHF_09*/"
                for var_name, var_type in self.outputs
            ],
            [  # b value check
                f"{var_type.test}({var_name}_a2d{var_type.get_bvalue}, {var_name}b_out);  /*UNQ_T2F_TFHF_10*/"
                for var_name, var_type in self.inputs
                if var_name in non_constant_inputs
            ],
            *[  # p value check
                [
                    f"{var_type.test}({var_name}_a2d{var_type.get_pvalue.format(i)}, {var_name}p{i}_out)"
                    f";  /*UNQ_T2F_TFHF_11*/"
                    for i in range(self.N)
                ]
                for var_name, var_type in self.outputs
            ]
        ))
        return self._test_f_(test_variant_name, "hforward", declarations_list,
                             derivative_loop, evaluations_list, comparisons_list)

    def _test_f_h_reverse_(self, test_variant_name: str, non_constant_inputs: non_constant_inputs_type) -> str:
        assert bool(non_constant_inputs)  # this might need to change later
        non_constant_input_types = {var_type for var_name, var_type in self.inputs if var_name in non_constant_inputs}
        declarations_list = (
                [  # initialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(itertools.chain(
                        # standard inputs
                        self._all_declaration_format_(
                            "{var_name}({var_name}_data)", self.inputs, var_type
                        ),
                        # "p" inputs (for the non-constant inputs)
                        *(self._non_constant_declaration_format_(
                            "{var_name}p{i}({var_name}p{i}_data)", self.inputs, var_type, non_constant_inputs, i=i
                        ) for i in range(self.N)),
                        # "b" outputs (FYI: if any of the inputs are non-constant, then all outputs are non-constant)
                        (self._all_declaration_format_(
                            "{var_name}b({var_name}b_data)", self.outputs, var_type
                        ) * bool(non_constant_inputs)),
                        # "h" outputs (FYI: if any of the inputs are non-constant, then all outputs are non-constant)
                        *((self._all_declaration_format_(
                            "{var_name}h{i}({var_name}h{i}_data)", self.outputs, var_type, i=i
                        ) * bool(non_constant_inputs))
                          for i in range(self.N))
                    )) +
                    ";  /*UNQ_T2F_TFHR_01*/"
                    for var_type in self.datatypes
                ] +
                [  # uninitialized declarations
                    f"{var_type.parent.alias} " +
                    ', '.join(itertools.chain(
                        # standard outputs
                        self._all_declaration_format_(
                            "{var_name}", self.outputs, var_type
                        ),
                        # "b" inputs (for the non-constant inputs)
                        self._non_constant_declaration_format_(
                            "{var_name}b", self.inputs, var_type, non_constant_inputs
                        )
                    )) +
                    ";  /*UNQ_T2F_TFHR_02*/"
                    # for var_type in self.output_types
                    # if var_type.requires_parent_declaration
                    for var_type in self.datatypes
                    if (var_type.requires_parent_declaration and
                        (var_type in non_constant_input_types or
                         var_type in self.output_types))
                ] +
                [  # dependent declarations (A.K.A. A2D initializations)
                    f"{var_type.alias} " +
                    ', '.join(
                        self._non_constant_declaration_format_(
                            "{var_name}_a2d({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_a2d({var_name}, {default_bvalue})",
                            self.inputs, var_type, non_constant_inputs
                        ) +
                        (self._all_declaration_format_(
                            "{var_name}_a2d({var_name}, {var_name}b)" if var_type.requires_parent_declaration else
                            "{var_name}_a2d({default_value}, {var_name}b)",
                            self.outputs, var_type
                        ) * bool(non_constant_inputs))
                    ) +
                    ";  /*UNQ_T2F_TFHR_03*/"
                    for var_type in self.datatypes
                    # Make sure this declaration isn't empty
                    if var_type in self.output_types or var_type in non_constant_input_types
                ]
        )

        derivative_loops_list = [
            [  # p value assignments
                f"{var_name}_a2d{var_type.get_pvalue.format(i)}{{rank{rank}}} = {var_name}p{i}{{rank{rank}}}"
                f";  /*UNQ_T2F_TFHR_04*/"
                for i in range(self.N)
                for var_name, var_type in self.inputs
                if (var_name in non_constant_inputs and
                    var_type.rank == rank)
            ] +
            [  # h value assignments
                f"{var_name}_a2d{var_type.get_hvalue.format(i)}{{rank{rank}}} = {var_name}h{i}{{rank{rank}}}"
                f";  /*UNQ_T2F_TFHR_05*/"
                for i in range(self.N)
                for var_name, var_type in self.outputs
                if var_type.rank == rank
            ]
            for rank in range(max(inp.rank for inp in non_constant_input_types.union(self.output_types)) + 1)
        ]
        derivative_loop = self._construct_derivative_loop_(derivative_loops_list)

        evaluations_list = [
            f'auto expr = A2D::{self.operation_name}({self._evaluation_signature_(non_constant_inputs)})'
            f';  /*UNQ_T2F_TFHR_06*/',
            f'expr.reverse();  /*UNQ_T2F_TFHR_07*/',
            f'expr.hforward();  /*UNQ_T2F_TFHR_08*/',
            f'expr.hreverse();  /*UNQ_T2F_TFHR_09*/',
        ]

        comparisons_list = list(itertools.chain(  # This depends on the type of test (e.g. passive, forward, etc...)
            [  # value check
                f"{var_type.test}({var_name}_a2d{var_type.get_value}, {var_name}_out);  /*UNQ_T2F_TFHR_10*/"
                for var_name, var_type in self.outputs
            ],
            [  # b value check
                f"{var_type.test}({var_name}_a2d{var_type.get_bvalue}, {var_name}b_out);  /*UNQ_T2F_TFHR_11*/"
                for var_name, var_type in self.inputs
                if var_name in non_constant_inputs
            ],
            *[  # p value check
                [
                    f"{var_type.test}({var_name}_a2d{var_type.get_pvalue.format(i)}, {var_name}p{i}_out)"
                    f";  /*UNQ_T2F_TFHR_12*/"
                    for i in range(self.N)
                ]
                for var_name, var_type in self.outputs
            ],
            *[  # h value check
                [
                    f"{var_type.test}({var_name}_a2d{var_type.get_hvalue.format(i)}, {var_name}h{i}_out, 1e-8)"
                    f";  /*UNQ_T2F_TFHR_13*/"
                    for i in range(self.N)
                ]
                for var_name, var_type in self.inputs
                if var_name in non_constant_inputs
            ]
        ))
        return self._test_f_(test_variant_name, "hreverse", declarations_list,
                             derivative_loop, evaluations_list, comparisons_list)

    def full_test_variant(self, input_data: input_data_type,
                          non_constant_inputs: non_constant_inputs_type) -> str:
        test_variant_name = self._test_variant_name_(non_constant_inputs)
        separator_string = '\n    '

        # generate outputs
        if non_constant_inputs:
            nc_inputs = [(var_name, var_type)
                         for var_name, var_type in self.inputs
                         if var_name in non_constant_inputs]
            # Include derivatives
            results = (
                    [
                        # calculated nominal values
                        var_type.initialize_assignment_format.format(var_name=f"{var_name}_out", *value)
                        for (var_name, var_type), value in zip(self.outputs, self.evaluate_at(input_data))
                    ] +
                    [
                        # calculated input derivatives (b values)
                        var_type.initialize_assignment_format.format(var_name=f"{var_name}b_out", *value)
                        for (var_name, var_type), value in
                        zip(nc_inputs, self.evaluate_input_derivatives_at(input_data, non_constant_inputs))
                    ] +
                    [
                        # calculated output derivatives (p values)
                        var_type.initialize_assignment_format.format(var_name=f"{var_name}p{i}_out", *value)
                        for i, output_derivatives in
                        enumerate(self.evaluate_output_derivatives_at(input_data, non_constant_inputs))
                        for (var_name, var_type), value in
                        zip(self.outputs, output_derivatives)
                    ] +
                    [
                        # calculated input second derivatives (h values)
                        var_type.initialize_assignment_format.format(var_name=f"{var_name}h{i}_out", *value)
                        for i, hessian_derivatives in
                        enumerate(self.evaluate_hessian_derivatives_at(input_data, non_constant_inputs))
                        for var_name, var_type, value in
                        hessian_derivatives
                    ]
            )
        else:
            # No derivatives necessary
            results = (
                [
                    # calculated nominal values
                    var_type.initialize_assignment_format.format(var_name=f"{var_name}_out", *value)
                    for (var_name, var_type), value in zip(self.outputs, self.evaluate_at(input_data))
                ]
            )

        test_string = f"""
class {test_variant_name} : public {self.operation_name} {{
protected:
    // Results
    {separator_string.join(results)}
}};  /*UNQ_T2F_FTV_01*/
{self._test_f_passive_(test_variant_name, non_constant_inputs)}"""

        if non_constant_inputs:
            return (test_string
                    + self._test_f_reverse_(test_variant_name, non_constant_inputs)
                    + self._test_f_h_forward_(test_variant_name, non_constant_inputs)
                    + self._test_f_h_reverse_(test_variant_name, non_constant_inputs)
                    )
        return test_string

    pass


class TestConstructor:
    """\
TestConstructor: An object that constructs a full test suite for functions defined by TestADFunction and/or \
TestA2DFunction objects.  Use the TestConstructor.construct method to write the test file.  The test can either \
be manually added to the test CMakeLists.txt file or an optional argument ("add_to_cmake") can be specified to \
automatically update the CMakeLists.txt file.

name: the name of the test file to be constructed.  The convention is to include all tests for objects in a \
"<name>.h" file in a single "test_<name>.cpp" file of the same name, however it may be advisable to split some header \
files between multiple test files.

packages_to_test: packages to include in addition to those automatically included for all test files.  Automatically \
included packages are "gtest/gtest.h", "a2dobjs.h", "a2dtypes.h", and "test_commons.h"

var_types: a list of all VarType, ADVarType, and A2DVarType objects directly or indirectly used by the \
TestADFunction and/or TestA2DFunction objects.

test_functions: a list of TestADFunction and/or TestA2DFunction objects describing the overloaded operations \
for which to write tests.
"""

    # Uses "/*UNQ_TC_<method acronym>_##*/" as the unique identifier comment for debugging purposes.

    name: str
    packages_to_test: List[str]
    var_types: List[Union[VarType, ADVarType]]
    test_functions: List[Union[TestADFunction, TestA2DFunction]]

    def __init__(self, name: str, packages_to_test: List[str], var_types: List[Union[VarType, ADVarType]],
                 test_functions: List[Union[TestADFunction, TestA2DFunction]]):
        self.name = name
        self.packages_to_test = packages_to_test
        self.var_types = var_types
        self.test_functions = test_functions
        pass

    @property
    def typedefs(self) -> Set[Union[VarType, ADVarType]]:
        # get var_types from test_functions and combine with given var_types
        return set(self.var_types).union(
            *[
                test_function.datatypes
                for test_function in self.test_functions
            ],
            *[
                {datatype.parent for datatype in test_function.datatypes}
                for test_function in self.test_functions
            ]
        )

    @staticmethod
    def __input_name_check__(var_name: str, all_io: dict, suffix: str):
        suffix_len = len(suffix)
        # make sure <var_name><suffix> is not in all_io
        if f"{var_name}{suffix}" in all_io:
            raise NamingConventionConflict(f'Illegal variable name "{var_name}{suffix}" when a \n'
                                           f'variable "{var_name}" is also defined.')
        # make sure <var_name>[:-suffix_len]## is not in all_io (when var_name ends with <suffix>##)
        numeric_stripped_var_name = var_name.rstrip('0123456789')
        if numeric_stripped_var_name[-suffix_len:] == suffix and numeric_stripped_var_name[:-suffix_len] in all_io:
            raise NamingConventionConflict(f'Illegal variable name "{var_name}" when a \n'
                                           f'variable "{numeric_stripped_var_name[:-suffix_len]}" is also defined.')
        pass

    def input_data(self) -> input_data_type:
        # gather all inputs and outputs and their types (and make sure there aren't conflicts)
        all_io = {}
        for test_function in self.test_functions:
            for var_name, var_type in itertools.chain(test_function.inputs, test_function.outputs):
                if var_name in all_io:
                    if all_io[var_name] != var_type:
                        raise VariableOverloadError(f'{var_name} assigned to multiple, non-equivalent types: '
                                                    f'{all_io[var_name]} and {var_type}')
                    self.__input_name_check__(var_name, all_io, 'b')
                    self.__input_name_check__(var_name, all_io, 'p')
                    self.__input_name_check__(var_name, all_io, 'h')
                all_io[var_name] = var_type
                pass
            pass

        a2d_inputs = {}
        a2d_outputs = {}
        for test_function in filter(TestA2DFunction.__instancecheck__, self.test_functions):
            a2d_inputs.update(
                {var_name: (var_type, max(test_function.N, a2d_inputs.get(var_name, (None, 0))[1]))
                 for var_name, var_type in test_function.inputs}
            )
            a2d_outputs.update(
                {var_name: (var_type, max(test_function.N, a2d_outputs.get(var_name, (None, 0))[1]))
                 for var_name, var_type in test_function.outputs}
            )
            pass

        all_inputs = {}
        for test_function in self.test_functions:
            all_inputs.update(test_function.inputs)
            pass

        # generate data
        data = {
            var_name: (var_type.parent.generator(), var_type)
            for var_name, var_type in all_inputs.items()
        }
        # generate b derivative data
        data.update({
            f"{var_name}b": (var_type.parent.generator(), var_type)
            for var_name, var_type in all_io.items()
        })
        # generate p derivative data for inputs
        data.update({
            f"{var_name}p{i}": (var_type.parent.generator(), var_type)
            for var_name, (var_type, N) in a2d_inputs.items()
            for i in range(N)
        })
        # generate h derivative data for outputs
        data.update({
            f"{var_name}h{i}": (var_type.parent.generator(), var_type)
            for var_name, (var_type, N) in a2d_outputs.items()
            for i in range(N)
        })
        return data  # {<var_name>: (<value>, <var_type>), ...}

    def __test_function__(self, input_data: input_data_type, test_function: TestADFunction):
        separator_string = '\n    '
        overloaded_data = []  # note: must combine overloaded_data with input_data to pass if implementing this
        return f"""
class {test_function.operation_name} : public {self.name}Test {{
protected:
    {separator_string.join(overloaded_data)}
}};  /*UNQ_TC_TF_01*/
{test_function.tests(input_data)}"""

    # Component Construction Methods:

    def _title_comment_(self, max_len=80):
        separator_string = '\n    '
        end_sep = (', and '
                   if len(self.packages_to_test) > 2 else
                   ' and '
                   if len(self.packages_to_test) == 2 else
                   '')
        name_string = ', '.join(self.packages_to_test[:-1]) + end_sep + self.packages_to_test[-1]
        base_string = (
            f"This is a set of automatically generated unit tests for {name_string} using Google Test framework.  "
            f"These tests were written on {datetime.date.today().isoformat()} using the {__package_name__} package."
        ).split(' ')
        base_string_lengths = list(map(len, base_string))
        split_indices = [0]
        for i in range(len(base_string)):
            if sum(base_string_lengths[split_indices[-1]: i]) + (i - split_indices[-1]) > max_len:
                split_indices.append(i - 1)
            pass
        split_indices.append(len(base_string))
        comment_string = separator_string.join(' '.join(base_string[i:j])
                                               for i, j in zip(split_indices, split_indices[1:]))
        return f"""\
/*
    {comment_string}
*/
"""

    def _include_statements_(self):
        includes = '\n'.join(f'#include "{pkg}"' for pkg in self.packages_to_test)
        return f"""
#include <gtest/gtest.h>

#include "a2dobjs.h"
#include "a2dtypes.h"
{includes}
#include "test_commons.h"
"""

    def _type_definitions_(self):
        # put typedefs in order based on precedence
        aliases = []
        definitions = []
        for typedef in self.typedefs:
            if any(map(typedef.alias_pattern.search, definitions)):
                position = list(map(bool, map(typedef.alias_pattern.search, definitions))).index(True)
                aliases.insert(position, typedef.alias)
                definitions.insert(position, typedef.definition)
            else:
                aliases.append(typedef.alias)
                definitions.append(typedef.definition)
            pass
        # format type definitions
        defs = '\n'.join(f'using {alias} = {definition};  /*UNQ_TC_TD_01*/'
                         for alias, definition in zip(aliases, definitions))
        return f"""
{defs}
"""

    def _testing_input_class_(self, input_data: input_data_type):
        separator_string = '\n    '
        data = [
            var_type.initialize_assignment_format.format(var_name=f"{var_name}_data", *value)
            for var_name, (value, var_type) in input_data.items()
        ]
        return f"""
class {self.name}Test : public ::testing::Test {{
protected:
    // Input Options:
    {separator_string.join(data)}
}};  /*UNQ_TC_TIC_01*/
"""

    def _all_tests_(self, input_data: input_data_type):
        return ''.join(
            self.__test_function__(input_data, test_function)
            for test_function in self.test_functions
        )

    def construct(self, destination: str = None, name_override: str = None, add_to_cmake: bool = False):
        input_data = self.input_data()
        file_str = (
                self._title_comment_() +
                self._include_statements_() +
                self._type_definitions_() +
                self._testing_input_class_(input_data) +
                self._all_tests_(input_data)
        )
        if destination is None:
            return file_str
        cpp_pattern = re.compile(r'\.cpp$')
        test_pattern = re.compile(r'^test_')
        name = test_pattern.sub('', cpp_pattern.sub('', (
            self.name if (name_override is None) else name_override).strip()))
        filename = f"test_{name}.cpp"
        separator = os.altsep if os.altsep in destination else os.sep
        full_filename = destination.rstrip(separator) + separator + filename
        with open(full_filename, 'w') as f:
            f.write(file_str)
            f.close()

        print(f'\n'
              f'{filename} successfully written to:\n'
              f'{destination}')
        if add_to_cmake:
            print("CMakeLists:")
            # dynamically add file to CMakeLists.txt
            if 'CMakeLists.txt' in os.listdir(destination):
                lines = [
                    f"add_executable({cpp_pattern.sub('', filename)} {filename})",
                    f"target_link_libraries({cpp_pattern.sub('', filename)} gtest_main)",
                    f"gtest_discover_tests({cpp_pattern.sub('', filename)})",
                ]
                with open(destination.rstrip(separator) + separator + 'CMakeLists.txt', 'r') as f:
                    initial_lines = [line.strip() for line in f.readlines()]
                    f.close()
                    pass

                command_not_found = False
                modified_lines = initial_lines[:]
                for new_line in lines:
                    command = new_line[:new_line.index('(') + 1]
                    if any(command in line for line in initial_lines):
                        if new_line not in initial_lines:
                            new_line_last_index = max(i for i, line in enumerate(modified_lines)
                                                      if command in line)
                            modified_lines.insert(new_line_last_index + 1, new_line)
                            pass
                    else:
                        w = UserWarning(
                            f'\n'
                            f'No commands matching "{command[:-1]}" currently in CMakeLists.txt file; it might be\n'
                            f'incompatible.  If you are certain the "{command[:-1]}" command can be added, please\n'
                            f'do so manually.')
                        warnings.warn(w)
                        command_not_found = True
                    pass

                if modified_lines != initial_lines:
                    with open(destination.rstrip(separator) + separator + 'CMakeLists.txt', 'w') as f:
                        f.write('\n'.join(modified_lines))
                        f.close()
                        print(f'CMakeLists.txt updated with {filename}.')
                        pass
                    pass
                elif command_not_found:
                    print(f'{filename} already has commands for all recognized commands\n'
                          f'in CMakeLists.txt, no update performed.')
                else:
                    print(f'{filename} already in CMakeLists.txt, no update required.')
                pass
            else:
                raise FileNotFoundError(f'Attempted to dynamically update "CMakeLists.txt", but no such file\n'
                                        f'was found in the given destination folder.')
            pass
        else:
            print(f'IMPORTANT: you must manually add this test suite to the test CMakeLists.txt file.')
            pass
        pass

    pass


if __name__ == '__main__':
    """Example construction of a2dvecops3d AD class tests"""


    def symmat():
        non_sym = numpy.random.random((3, 3))
        return (non_sym + non_sym.T).flatten()


    i_t = VarType(alias='I',
                  shortname='Int',
                  definition='int',
                  test='expect_val_eq',
                  initialize_assignment_format='const T {var_name} = {0};',
                  rank=0,
                  generator=lambda: numpy.random.randint(-256, 256, 1))
    t_t = VarType(alias='T',
                  shortname='Scalar',
                  definition='double',
                  test='expect_val_eq',
                  initialize_assignment_format='const T {var_name} = {0};',
                  rank=0,
                  generator=lambda: numpy.random.random(1))
    vec_t = VarType(alias='Vec_t',
                    shortname='Vec',
                    definition='A2D::Vec<T, 3>',
                    test='expect_vec_eq<3>',
                    initialize_assignment_format='const T {var_name}[3] = {{{0}, {1}, {2}}};',
                    rank=1,
                    generator=lambda: numpy.random.random(3))
    mat_t = VarType(alias='Mat_t',
                    shortname='Mat',
                    definition='A2D::Mat<T, 3, 3>',
                    test='expect_mat_eq<3, 3>',
                    initialize_assignment_format=
                    'const T {var_name}[9] = {{{0}, {1}, {2}, \n{3}, {4}, {5}, \n{6}, {7}, {8}}};',
                    rank=2,
                    generator=symmat)
    advec_t = ADVarType(alias='ADVec_t',
                        shortname='ADVec',
                        definition='A2D::ADVec<Vec_t>',
                        get_value='.value()',
                        get_bvalue='.bvalue()',
                        parent=vec_t)
    adscalar_t = ADVarType(alias='ADScalar_t',
                           shortname='ADScalar',
                           definition='A2D::ADScalar<T>',
                           get_value='.value',
                           get_bvalue='.bvalue',
                           parent=t_t,
                           initialization_defaults=(0, 0))
    admat_t = ADVarType(alias='ADMat_t',
                        shortname='ADMat',
                        definition='A2D::ADMat<Mat_t>',
                        get_value='.value()',
                        get_bvalue='.bvalue()',
                        parent=mat_t)

    tfs = [
        TestADFunction("Vec3Norm",
                       inputs=[('x', advec_t), ],
                       outputs=[('a', adscalar_t), ],
                       operation=lambda x: (numpy.array([numpy.dot(x, x) ** 0.5]),),
                       variants=[{'x'}, set()]),
        TestADFunction("Vec3Scale",
                       inputs=[('x', advec_t),
                               ('a', adscalar_t), ],
                       outputs=[('v', advec_t), ],
                       operation=lambda x, a: (numpy.multiply(a, x),),
                       variants=[{'x', 'a'}, {'x'}, {'a'}, set()]),
        TestADFunction("Vec3Dot",
                       inputs=[('x', advec_t),
                               ('y', advec_t), ],
                       outputs=[('a', adscalar_t), ],
                       operation=lambda x, y: (numpy.array([numpy.dot(x, y)]),),
                       variants=[{'x', 'y'}, {'x'}, set()]),
        TestADFunction("Vec3Normalize",
                       inputs=[('x', advec_t), ],
                       outputs=[('v', advec_t), ],
                       operation=lambda x: (numpy.multiply(1 / (numpy.dot(x, x) ** 0.5), x),),
                       variants=[{'x'}, set()]),
        TestADFunction("Vec3ScaleSymmetricOuterProduct",
                       inputs=[('a', adscalar_t),
                               ('x', advec_t), ],
                       outputs=[('S', admat_t), ],
                       operation=lambda a, x: (numpy.multiply(a, numpy.outer(x, x)).flatten(),),
                       variants=[{'a', 'x'}, {'a'}, {'x'}, set()]),
    ]

    tc = TestConstructor(name="vecops3d_jg",
                         packages_to_test=["a2dvecops3d.h"],
                         var_types=[t_t, i_t, vec_t, mat_t, advec_t, adscalar_t, admat_t],
                         test_functions=tfs)

    tc.construct(destination=os.getcwd(),
                 add_to_cmake=True)
    pass

if __name__ == '__main__':
    """Example construction of a2dvecops3d A2D class tests."""


    def symmat():
        non_sym = numpy.random.random((3, 3))
        return (non_sym + non_sym.T).flatten()


    N_ = 4
    t_t = VarType(alias='T',
                  shortname='Scalar',
                  definition='double',
                  test='expect_val_eq',
                  initialize_assignment_format='const T {var_name} = {0};',
                  rank=0,
                  generator=lambda: numpy.random.random(1))
    vec_t = VarType(alias='Vec_t',
                    shortname='Vec',
                    definition='A2D::Vec<T, 3>',
                    test='expect_vec_eq<3>',
                    initialize_assignment_format='const T {var_name}[3] = {{{0}, {1}, {2}}};',
                    rank=1,
                    generator=lambda: numpy.random.random(3))
    mat_t = VarType(alias='Mat_t',
                    shortname='Mat',
                    definition='A2D::Mat<T, 3, 3>',
                    test='expect_mat_eq<3, 3>',
                    initialize_assignment_format=
                    'const T {var_name}[9] = {{{0}, {1}, {2}, \n{3}, {4}, {5}, \n{6}, {7}, {8}}};',
                    rank=2,
                    generator=symmat)
    a2dvec_t = A2DVarType(alias='A2DVec_t',
                          shortname='A2DVec',
                          definition='A2D::A2DVec<{N}, Vec_t>',
                          get_value='.value()',
                          get_bvalue='.bvalue()',
                          get_pvalue='.pvalue({0})',
                          get_hvalue='.hvalue({0})',
                          parent=vec_t,
                          N=N_)
    a2dscalar_t = A2DVarType(alias='A2DScalar_t',
                             shortname='A2DScalar',
                             definition='A2D::A2DScalar<{N}, T>',
                             get_value='.value',
                             get_bvalue='.bvalue',
                             get_pvalue='.pvalue[{0}]',
                             get_hvalue='.hvalue[{0}]',
                             parent=t_t,
                             N=N_,
                             initialization_defaults=(0, 0))
    a2dmat_t = A2DVarType(alias='A2DMat_t',
                          shortname='A2DMat',
                          definition='A2D::A2DMat<{N}, Mat_t>',
                          get_value='.value()',
                          get_bvalue='.bvalue()',
                          get_pvalue='.pvalue({0})',
                          get_hvalue='.hvalue({0})',
                          parent=mat_t,
                          N=N_)

    tfs = [
        TestA2DFunction("Vec3Norm",
                        inputs=[('x', a2dvec_t), ],
                        outputs=[('a', a2dscalar_t), ],
                        operation=lambda x: (numpy.array([numpy.dot(x, x) ** 0.5]),),
                        variants=[{'x'}]),
        TestA2DFunction("Vec3Scale",
                        inputs=[('x', a2dvec_t),
                                ('a', a2dscalar_t), ],
                        outputs=[('v', a2dvec_t), ],
                        operation=lambda x, a: (numpy.multiply(a, x),),
                        variants=[{'x', 'a'}, {'x'}, {'a'}]),
        TestA2DFunction("Vec3Axpy",
                        inputs=[('a', a2dscalar_t),
                                ('x', a2dvec_t),
                                ('y', a2dvec_t), ],
                        outputs=[('v', a2dvec_t), ],
                        operation=lambda a, x, y: (numpy.multiply(a, x) + y,),
                        variants=[{'a', 'x', 'y'}, {'a', 'x'}, {'x', 'y'}, {'a', 'y'}, {'a'}, {'x'}, {'y'}]),
        TestA2DFunction("Vec3Dot",
                        inputs=[('x', a2dvec_t),
                                ('y', a2dvec_t), ],
                        outputs=[('a', a2dscalar_t), ],
                        operation=lambda x, y: (numpy.array([numpy.dot(x, y)]),),
                        variants=[{'x', 'y'}, {'x'}, {'y'}]),
        TestA2DFunction("Vec3Cross",
                        inputs=[('x', a2dvec_t),
                                ('y', a2dvec_t), ],
                        outputs=[('v', a2dvec_t), ],
                        operation=lambda x, y: (numpy.cross(x, y),),
                        variants=[{'x', 'y'}, {'x'}, {'y'}]),
        TestA2DFunction("Vec3Normalize",
                        inputs=[('x', a2dvec_t), ],
                        outputs=[('v', a2dvec_t), ],
                        operation=lambda x: (numpy.multiply(1 / (numpy.dot(x, x) ** 0.5), x),),
                        variants=[{'x'}]),
        TestA2DFunction("Vec3ScaleSymmetricOuterProduct",
                        inputs=[('a', a2dscalar_t),
                                ('x', a2dvec_t), ],
                        outputs=[('S', a2dmat_t), ],
                        operation=lambda a, x: (numpy.multiply(a, numpy.outer(x, x)).flatten(),),
                        variants=[{'a', 'x'}, {'a'}, {'x'}]),
    ]

    tc = TestConstructor(name="vecops3d_a2d",
                         packages_to_test=["a2dvecops3d.h"],
                         var_types=[t_t, vec_t, mat_t, a2dvec_t, a2dscalar_t, a2dmat_t],
                         test_functions=tfs)

    tc.construct(destination=os.getcwd(),
                 add_to_cmake=True)
    pass
