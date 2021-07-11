from amuse.community import CodeInterface, LiteratureReferencesMixIn
from amuse.community.interface.gd import GravitationalDynamics, GravitationalDynamicsInterface
from amuse.rfi.core import legacy_function, LegacyFunctionSpecification
from amuse.units import nbody_system

class FalconInterface(CodeInterface, GravitationalDynamicsInterface, LiteratureReferencesMixIn):
    '''
    Fast-multipole Poisson solver "falcON"

    .. [#] Dehnen W., 2000, ApJL, 536, L39
    .. [#] Dehnen W., 2002, J.Comp.Phys., 179, 27
    '''
    include_headers = ['worker_code.h']

    def __init__ (self, **kwargs):
        CodeInterface.__init__(self, name_of_the_worker="falcon_worker", **kwargs),
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def set_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='float64', direction=function.IN, description = 'timestep')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_epsilon():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='float64', direction=function.OUT,
            description = 'softening length', unit = nbody_system.length)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_epsilon():
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon', dtype='float64', direction=function.IN,
            description = 'softening length', unit = nbody_system.length)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_individual_epsilon():
        function = LegacyFunctionSpecification()
        function.addParameter('flag', dtype='bool', direction=function.IN,
            description = 'enable/disable individual softening lengths')
        function.result_type = 'int32'
        return function


class Falcon(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **kwargs):
        GravitationalDynamics.__init__(self, FalconInterface(**kwargs), convert_nbody, **kwargs)
        for name, value in kwargs.items():
            print('setting %s to %s' % (name, value))
            setattr(self.parameters, name, value)

    def define_parameters(self, handler):
        handler.add_method_parameter(
            None,
            'set_individual_epsilon',
            'individual_epsilon',
            'use individual softening lengths for each particle (their radii)',
        )
        handler.add_method_parameter(
            'get_epsilon',
            'set_epsilon',
            'epsilon',
            'softening length',
        )
        handler.add_method_parameter(
            'get_eps2',
            'set_eps2',
            'epsilon_squared',
            'squared softening length',
        )
        handler.add_method_parameter(
            'get_time_step',
            'set_time_step',
            'timestep',
            'constant timestep for iteration',
        )

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        handler.add_method(
            "get_time_step",
            (),
            (nbody_system.time, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_time_step",
            (nbody_system.time, ),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_eps2",
            (),
            (nbody_system.length**2, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_eps2",
            (nbody_system.length**2, ),
            (handler.ERROR_CODE,)
        )
