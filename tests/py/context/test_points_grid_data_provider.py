import pytest
from pydantic import ValidationError
from mat3ra.esse.models.context_providers_directory.points_grid_data_provider import (
    GridMetricType,
    PointsGridDataProviderSchema,
)
from mat3ra.made.material import Material
from mat3ra.wode.context.providers import PointsGridDataProvider
from mat3ra.wode.context.providers.points_grid_data_provider import DEFAULT_KPPRA


def _material_stub(number_of_atoms, reciprocal_vector_ratios):
    a, b, c = (2.0 / ratio for ratio in reciprocal_vector_ratios)
    return Material.create(
        {
            "name": "test",
            "lattice": {
                "a": a,
                "b": b,
                "c": c,
                "alpha": 90,
                "beta": 90,
                "gamma": 90,
                "type": "ORC",
                "units": {"length": "angstrom", "angle": "degree"},
            },
            "basis": {
                "elements": [{"id": i, "value": "Si"} for i in range(number_of_atoms)],
                "coordinates": [{"id": i, "value": [i * 0.1] * 3} for i in range(number_of_atoms)],
                "units": "crystal",
            },
        }
    )

# Test data constants
DIMENSIONS_DEFAULT = [1, 1, 1]
DIMENSIONS_CUSTOM = [1, 2, 3]
SHIFTS_DEFAULT = [0.0, 0.0, 0.0]
SHIFTS_CUSTOM = [0.5, 0.5, 0.5]
DIVISOR_DEFAULT = 1
DIVISOR_CUSTOM = 2
GRID_METRIC_TYPE_DEFAULT = GridMetricType.KPPRA
NUMBER_OF_ATOMS_DEFAULT = 1
RATIOS_DEFAULT = [1.0, 1.0, 1.0]
GRID_METRIC_VALUE_DERIVED = (
    DIMENSIONS_CUSTOM[0] * DIMENSIONS_CUSTOM[1] * DIMENSIONS_CUSTOM[2] * NUMBER_OF_ATOMS_DEFAULT
)

# Expected data structures
KGRID_DATA = {
    "kgrid": {
        "dimensions": DIMENSIONS_CUSTOM,
        "shifts": SHIFTS_DEFAULT,
        "gridMetricType": GRID_METRIC_TYPE_DEFAULT,
        "gridMetricValue": GRID_METRIC_VALUE_DERIVED,
        "preferGridMetric": False,
        "reciprocalVectorRatios": RATIOS_DEFAULT,
    },
    "isKgridEdited": True,
}

KGRID_TEMPLATE_DATA = {
    "kgrid": {
        "dimensions": ["{{N_k}}", "{{N_k}}", "{{N_k}}"],
        "shifts": SHIFTS_DEFAULT,
        "gridMetricType": GRID_METRIC_TYPE_DEFAULT,
        "preferGridMetric": False,
        "gridMetricValue": DEFAULT_KPPRA,
        "reciprocalVectorRatios": [1.0, 0.667, 0.5],
    },
    "isKgridEdited": True,
    "isUsingJinjaVariables": True,
}


@pytest.mark.parametrize(
    "init_params,expected_dimensions,expected_shifts,expected_divisor",
    [
        (
            {"dimensions": DIMENSIONS_CUSTOM, "material": _material_stub(NUMBER_OF_ATOMS_DEFAULT, RATIOS_DEFAULT)},
            DIMENSIONS_CUSTOM,
            SHIFTS_DEFAULT,
            DIVISOR_DEFAULT,
        ),
    ],
)
def test_points_grid_data_provider_initialization(init_params, expected_dimensions, expected_shifts, expected_divisor):
    kgrid_context_provider_relax = PointsGridDataProvider(**init_params)

    assert kgrid_context_provider_relax.dimensions == expected_dimensions
    assert kgrid_context_provider_relax.shifts == expected_shifts
    assert kgrid_context_provider_relax.divisor == expected_divisor


@pytest.mark.parametrize(
    "init_params,expected_data",
    [
        (
            {"dimensions": DIMENSIONS_CUSTOM, "material": _material_stub(NUMBER_OF_ATOMS_DEFAULT, RATIOS_DEFAULT)},
            KGRID_DATA,
        ),
    ],
)
def test_points_grid_data_provider_get_data(init_params, expected_data):
    kgrid_context_provider = PointsGridDataProvider(**init_params)
    actual_data = kgrid_context_provider.get_data()
    assert actual_data == expected_data["kgrid"]


@pytest.mark.parametrize(
    "init_params,expected_data",
    [
        (
            {
                "dimensions": DIMENSIONS_CUSTOM,
                "is_edited": True,
                "material": _material_stub(NUMBER_OF_ATOMS_DEFAULT, RATIOS_DEFAULT),
            },
            KGRID_DATA,
        ),
    ],
)
def test_points_grid_data_provider_yield_data(init_params, expected_data):
    kgrid_context_provider = PointsGridDataProvider(**init_params)
    actual_data = kgrid_context_provider.yield_data()
    assert actual_data == expected_data


def test_default_data_conforms_to_esse_schema():
    """Emitted data must validate against ESSE and carry no extra keys -- guards field drift."""
    material = _material_stub(NUMBER_OF_ATOMS_DEFAULT, RATIOS_DEFAULT)
    data = PointsGridDataProvider(dimensions=DIMENSIONS_CUSTOM, material=material).get_data()

    PointsGridDataProviderSchema.model_validate(data)  # raises on missing/wrong required fields
    assert set(data).issubset(set(PointsGridDataProviderSchema.model_fields))  # no subclass-only keys


def test_points_grid_data_provider_get_reciprocal_vector_ratios_from_provider_data():
    provider = PointsGridDataProvider(reciprocal_vector_ratios=[1.0, 0.667, 0.5])

    assert provider.get_reciprocal_vector_ratios() == [1.0, 0.667, 0.5]


def test_points_grid_data_provider_yield_data_with_overrides():
    provider = PointsGridDataProvider()

    actual_data = provider.yield_data_with_overrides(
        dimensions=["{{N_k}}", "{{N_k}}", "{{N_k}}"],
        reciprocal_vector_ratios=[1.0, 0.667, 0.5],
        is_using_jinja_variables=True,
    )

    assert actual_data == KGRID_TEMPLATE_DATA


def test_points_grid_data_provider_get_reciprocal_vector_ratios_from_context():
    provider = PointsGridDataProvider(
        context={
            "kgrid": {"reciprocalVectorRatios": [1.0, 0.8, 0.6]},
            "isKgridEdited": True,
        }
    )

    assert provider.get_reciprocal_vector_ratios() == [1.0, 0.8, 0.6]


def test_points_grid_data_provider_raises_when_atom_count_is_unknown():
    with pytest.raises(ValidationError, match="KPPRA"):
        PointsGridDataProvider(dimensions=[4, 4, 4], isEdited=True)


def test_material_derived_ratios_land_in_the_schema_field():
    material = _material_stub(number_of_atoms=2, reciprocal_vector_ratios=[1.0, 0.5, 0.25])

    provider = PointsGridDataProvider(dimensions=[4, 4, 4], material=material, isEdited=True)

    assert provider.reciprocalVectorRatios == [1.0, 0.5, 0.25]
    schema_only = provider.model_dump(
        by_alias=True, exclude_none=True, include=set(PointsGridDataProviderSchema.model_fields)
    )
    assert schema_only["reciprocalVectorRatios"] == [1.0, 0.5, 0.25]


def test_points_grid_data_provider_derives_grid_metric_value_and_ratios_from_material():
    material = _material_stub(number_of_atoms=2, reciprocal_vector_ratios=[1.0, 0.5, 0.25])

    data = PointsGridDataProvider(dimensions=[4, 4, 4], material=material, isEdited=True).get_data()

    assert data["gridMetricValue"] == 4 * 4 * 4 * 2
    assert data["reciprocalVectorRatios"] == [1.0, 0.5, 0.25]


def test_points_grid_data_provider_respects_explicit_grid_metric_value():
    provider = PointsGridDataProvider(dimensions=[4, 4, 4], isEdited=True, gridMetricValue=999)

    assert provider.get_data()["gridMetricValue"] == 999


def test_points_grid_data_provider_rounds_ratios_to_three_figures_like_js():
    material = _material_stub(2, [1.0, 0.8164965809277261, 0.6123724356957945])

    provider = PointsGridDataProvider(dimensions=[4, 4, 4], material=material, isEdited=True)

    assert provider.get_data()["reciprocalVectorRatios"] == [1.0, 0.816, 0.612]


def test_points_grid_data_provider_rejects_a_non_material():
    with pytest.raises(ValidationError):
        PointsGridDataProvider(dimensions=[4, 4, 4], material={"basis": 1, "lattice": 2})


def test_points_grid_data_provider_untouched_default_keeps_sentinel():
    provider = PointsGridDataProvider()

    assert provider.get_data()["gridMetricValue"] == DEFAULT_KPPRA


def test_standata_material_drives_derivation():
    from mat3ra.standata.materials import Materials

    material = Material.create(Materials.get_by_name_first_match("Silicon"))

    data = PointsGridDataProvider(dimensions=[4, 4, 4], material=material, isEdited=True).get_data()

    assert data["gridMetricValue"] == 4 * 4 * 4 * material.basis.number_of_atoms
    assert len(data["reciprocalVectorRatios"]) == 3
