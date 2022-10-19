import datetime
import inspect
import itertools
import os
import warnings
from typing import List, Dict, Callable, Any, Set, Tuple, Union, FrozenSet
import numpy

__package_name__ = __package__ if __package__ else "A2DTestConstructor"
var_generator_type = Callable[[], numpy.ndarray]


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
    generator: var_generator_type

    def __init__(self,
                 alias: str,
                 shortname: str,
                 definition: str,
                 test: str,
                 initialize_assignment_format: str,
                 generator: var_generator_type,
                 ):
        self.alias = alias
        self.shortname = shortname
        self.definition = definition
        self.test = test
        self.initialize_assignment_format = initialize_assignment_format
        self.generator = generator
        pass

    def __new__(cls,
                alias: str,
                shortname: str,
                definition: str,
                test: str,
                initialize_assignment_format: str,
                generator: var_generator_type,
                ):
        if alias in cls.__instances_by_alias__:
            equivalence = [getattr(cls.__instances_by_alias__[alias], key) == locals()[key]
                           for key in cls.__comparison_keys__()]
            if all(equivalence):
                return cls.__instances_by_alias__[alias]
            raise ValueError(f'Two different VarType objects with the same alias are not allowed.  There is already\n'
                             f'a VarType object with the alias "{alias}" with different attributes than were defined.')
        if alias in ADVarType.__instances_by_alias__:
            raise ValueError(f'A VarType object cannot have the same alias as a ADVarType object.')

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

definition: the specific C or C++ type used to define the type as alias.  For example, "A2D::Vec<T, 3>" or "double".

get_value: the method call operation used to obtain the base value for the type.  Usually ".value" or ".value()".

get_bvalue: the method call operation used to obtain the seed or bar value for the type.  Usually ".bvalue" or \
".bvalue()".

parent: the VarType object this ADVarType object inherits from.  Both the value and bvalue of object of this type \
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
    # initialization_format: str
    parent: VarType
    # requires_parent_declaration: bool
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
                 # initialization_format: str,
                 parent: VarType,
                 # requires_parent_declaration: bool,
                 initialization_defaults: Union[Tuple[Any, Any], Tuple[()]] = (),
                 ):
        self.alias = alias
        self.shortname = shortname
        self.definition = definition
        self.get_value = get_value
        self.get_bvalue = get_bvalue
        # self.initialization_format = initialization_format
        self.parent = parent
        # self.requires_parent_declaration = requires_parent_declaration
        self.initialization_defaults = initialization_defaults
        pass

    def __new__(cls,
                alias: str,
                shortname: str,
                definition: str,
                get_value: str,
                get_bvalue: str,
                # initialization_format: str,
                parent: VarType,
                # requires_parent_declaration: bool,
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

        obj = object.__new__(cls)
        cls.__instances_by_alias__[alias] = obj
        return obj

    @property
    def test(self):
        return self.parent.test

    @property
    def initialize_assignment_format(self):
        return self.parent.initialize_assignment_format

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


class VariableOverloadError(Exception):
    """Multiple types (VarTypes and/or ADVarTypes) assigned to same variable name."""
    pass


class NamingConventionConflict(Exception):
    """Collision of variable names for set naming convention."""
    pass


class VariableNameOverloadError(Exception):
    """Multiple variables with the same name."""
    pass


io_type = List[Tuple[str, ADVarType]]  # [(<var_name>, <var_type>), ...]
non_constant_inputs_type = FrozenSet[str]  # frozenset(<var_name>, ...)  <-These are the active inputs
non_constant_inputs_input_type = Set[str]  # {<var_name>, ...}  <-These are the active inputs
variants_type = List[non_constant_inputs_type]
variants_input_type = Union[List[non_constant_inputs_input_type],
                            Set[non_constant_inputs_type],
                            variants_type]
input_data_type = Dict[str, Tuple[Any, ADVarType]]  # {<var_name>: (<value>, <var_type>), ...}
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
IMPORTANT: all outputs must be returned in a tuple of non-scalar numpy.ndarray objects, even scalars.  Scalar outputs \
should be contained in a 1-dimensional, length 1 numpy.ndarray.  For example, an operation that has a single scalar \
output should return "(numpy.array([scalar_output]), )".

variants: (optional) a list of sets of input variable names defining the variations of the function to test.  Each \
set of input variable names defines a test case where the given inputs are assumed to be differentiable, all other \
inputs are assumed to be constants.  If left undefined, then the power set of variable names from "inputs" is used.
"""

    # Uses "/*UNQ_T1F_<method acronym>_##*/" as the unique identifier comment for debugging purposes.

    __h__ = 1e-30

    operation_name: str
    operation: operation_type
    inputs: io_type  # [(<var_name>, <var_type>), ...]
    outputs: io_type  # [(<var_name>, <var_type>), ...]
    variants: variants_type

    def __init__(self,
                 operation_name: str,
                 inputs: io_type,
                 outputs: io_type,
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
                                            f'TestADFunction({self.operation_name}, ...):\n'
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

    def __inspect_operation__(self):
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

    def __expected_variants__(self):
        input_names = [var_name for var_name, var_type in self.inputs]
        input_name_power_set = set(
            map(frozenset, itertools.chain.from_iterable(itertools.combinations(input_names, r)
                                                         for r in range(len(input_names) + 1)))
        )
        return input_name_power_set

    def evaluate_at(self, input_data: input_data_type):
        kwargs = {var_name: input_data[var_name][0]
                  for var_name, var_type in self.inputs}
        return self.operation(**kwargs)

    def evaluate_output_derivatives_at(self, input_data: input_data_type,
                                       non_constant_inputs: non_constant_inputs_type):
        kwargs = {
            var_name: (
                input_data[var_name][0] + (1j * self.__h__ * input_data[f"{var_name}b"][0])
                if var_name in non_constant_inputs else
                input_data[var_name][0]
            )
            for var_name, var_type in self.inputs
        }
        return [numpy.imag(x) / self.__h__ for x in self.operation(**kwargs)]

    def evaluate_input_derivatives_at(self, input_data: input_data_type,
                                      non_constant_inputs: non_constant_inputs_type):
        output_names = next(
            zip(*self.outputs))
        dimension_vectors = {
            var_name: [numpy.reshape(dim, numpy.shape(input_data[var_name][0]))
                       for dim in numpy.identity(numpy.size(input_data[var_name][0]))]
            for var_name, var_type in self.inputs
            if var_name in non_constant_inputs
        }
        partial_derivatives = {
            candidate_var_name: [
                {output_name: numpy.imag(x) / self.__h__
                 for output_name, x in
                 zip(output_names, self.operation(**{
                     var_name: (
                         input_data[var_name][0] + (1j * self.__h__ * dimension_vector)
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
    def _all_declaration_format_(string_format: str, io_vars: io_type, var_type: ADVarType) -> List[str]:
        defaults = dict(zip(['default_value', 'default_bvalue'], var_type.initialization_defaults))
        return [string_format.format(var_name=var_name, **defaults)
                for var_name, io_var_type in io_vars
                if io_var_type == var_type]  # This might need modification

    @staticmethod
    def _non_constant_declaration_format_(string_format: str, io_vars: io_type, var_type: ADVarType,
                                          non_constant_inputs: non_constant_inputs_type) -> List[str]:
        defaults = dict(zip(['default_value', 'default_bvalue'], var_type.initialization_defaults))
        return [string_format.format(var_name=var_name, **defaults)
                for var_name, io_var_type in io_vars
                if (io_var_type == var_type and  # This might need modification
                    var_name in non_constant_inputs)]

    def _evaluation_signature_(self, non_constant_inputs: non_constant_inputs_type):
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


class TestConstructor:
    """\
TestConstructor: An object that constructs a full test suite for functions defined by TestADFunction objects.  Use the \
.construct method to write the test file.  The test must be manually added to the test CMakeLists.txt file.

name: the name of the test file to be constructed.  The convention is to include all tests for objects in a \
"<name>.h" file in a single "test_<name>.cpp" file of the same name.

packages: packages to include in addition to those automatically included for all test files.  Automatically included \
packages are "gtest/gtest.h", "a2dobjs.h", "a2dtypes.h", and "test_commons.h"

var_types: a list of all VarType and ADVarType objects directly or indirectly used by the TestADFunction objects.

test_functions: a list of TestADFunction objects describing the overloaded operations to be tested.
"""

    # Uses "/*UNQ_TC_<method acronym>_##*/" as the unique identifier comment for debugging purposes.

    name: str
    packages: List[str]
    var_types: List[Union[VarType, ADVarType]]
    test_functions: List[TestADFunction]

    def __init__(self, name: str, packages_to_test: List[str], var_types: List[Union[VarType, ADVarType]],
                 test_functions: List[TestADFunction]):
        self.name = name
        self.packages = packages_to_test
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
        # make sure <var_name>[:-suffix_len] is not in all_io (when var_name ends with <suffix>)
        if var_name[-suffix_len:] == suffix and var_name[:-suffix_len] in all_io:
            raise NamingConventionConflict(f'Illegal variable name "{var_name}" when a \n'
                                           f'variable "{var_name[:-suffix_len]}" is also defined.')
        pass

    def input_data(self) -> input_data_type:
        # gather all inputs and outputs and their types (and make sure there aren't conflicts)
        all_io = {}
        for test_function in self.test_functions:
            for var_name, var_type in itertools.chain(test_function.inputs, test_function.outputs):
                if var_name in all_io:
                    if all_io[var_name] != var_type:
                        raise VariableOverloadError(f'{var_name} assigned multiple non-equivalent types: '
                                                    f'{all_io[var_name]} and {var_type}')
                    self.__input_name_check__(var_name, all_io, 'b')
                    self.__input_name_check__(var_name, all_io, 'p')
                    self.__input_name_check__(var_name, all_io, 'h')
                all_io[var_name] = var_type
                pass
            pass

        # generate data
        data = {
            var_name: (var_type.parent.generator(), var_type)
            for var_name, var_type in all_io.items()
        }
        # generate b derivative data
        data.update({
            f"{var_name}b": (var_type.parent.generator(), var_type)
            for var_name, var_type in all_io.items()
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
                   if len(self.packages) > 2 else
                   ' and '
                   if len(self.packages) == 2 else
                   '')
        name_string = ', '.join(self.packages[:-1]) + end_sep + self.packages[-1]
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
        includes = '\n'.join(f'#include "{pkg}"' for pkg in self.packages)
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
            if any(typedef.alias in typename for typename in definitions):
                position = [typedef.alias in typename for typename in definitions].index(True)
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
            # map(self.__test_function__, self.test_functions)
            self.__test_function__(input_data, test_function)
            for test_function in self.test_functions
        )

    def construct(self, destination: str = None, name_override: str = None):
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
        name = (self.name if (name_override is None) else name_override
                ).strip().rstrip('.cpp').lstrip('test_')
        filename = f"test_{name}.cpp"
        separator = os.altsep if os.altsep in destination else os.sep
        full_filename = destination.rstrip(separator) + separator + filename
        with open(full_filename, 'w') as f:
            f.write(file_str)
            f.close()

        print(f'\n{filename} successfully written to:\n'
              f'{destination}\n'
              f'IMPORTANT: you must manually add this test suite to the test CMakeLists.txt file.')
        pass

    pass


if __name__ == '__main__':
    """Construction of vector_operation_development tests"""


    def symmat():
        non_sym = numpy.random.random((3, 3))
        return (non_sym + non_sym.T).flatten()


    t_t = VarType('T', 'Scalar', 'double', 'expect_val_eq', 'const T {var_name} = {0};', lambda: numpy.random.random(1))
    i_t = VarType('I', 'Int', 'int', 'expect_val_eq', 'const T {var_name} = {0};',
                  lambda: numpy.random.randint(-256, 256, 1))
    vec_t = VarType('Vec_t', 'Vec', 'A2D::Vec<T, 3>', 'expect_vec_eq<3>', 'const T {var_name}[3] = {{{0}, {1}, {2}}};',
                    lambda: numpy.random.random(3))
    mat_t = VarType('Mat_t', 'Mat', 'A2D::Mat<T, 3, 3>', 'expect_mat_eq<3, 3>',
                    'const T {var_name}[9] = {{{0}, {1}, {2}, \n{3}, {4}, {5}, \n{6}, {7}, {8}}};',
                    symmat,
                    )
    advec_t = ADVarType('ADVec_t', 'ADVec', 'A2D::ADVec<Vec_t>', '.value()', '.bvalue()', vec_t)
    adscalar_t = ADVarType('ADScalar_t', 'ADScalar', 'A2D::ADScalar<T>', '.value', '.bvalue', t_t, (0, 0))
    admat_t = ADVarType('ADMat_t', 'ADMat', 'A2D::ADMat<Mat_t>', '.value()', '.bvalue()', mat_t)

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

    tc = TestConstructor("vecops3d_jg", ["a2dvecops3d.h"],
                         [t_t, i_t, vec_t, advec_t, adscalar_t],
                         tfs)

    tc.construct(os.getcwd())
    pass
