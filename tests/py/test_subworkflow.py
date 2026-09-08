import pytest
from mat3ra.ade.application import Application
from mat3ra.made.material import Material
from mat3ra.mode.method import Method
from mat3ra.mode.model import Model
from mat3ra.mode.models.dft import DFTModel
from mat3ra.standata.applications import ApplicationStandata
from mat3ra.standata.materials import Materials
from mat3ra.standata.workflows import WorkflowStandata

from mat3ra.wode import Subworkflow, Unit, Workflow, ExecutionUnit
from mat3ra.wode.context.providers import PointsGridDataProvider


SUBWORKFLOW_NAME = "Total Energy"
SUBWORKFLOW_APPLICATION = Application(**ApplicationStandata.get_by_name_first_match("espresso"))
SUBWORKFLOW_METHOD = Method(type="pseudopotential", subtype="us")
SUBWORKFLOW_MODEL = Model(type="dft", subtype="gga", method=SUBWORKFLOW_METHOD)
SUBWORKFLOW_PROPERTIES = ["total_energy", "pressure"]
WORKFLOW_STANDATA = WorkflowStandata()
APPLICATION_ESPRESSO = ApplicationStandata.get_by_name_first_match("espresso")["name"]
DEFAULT_WF_NAME = WORKFLOW_STANDATA.get_default()["name"]

UNIT_CONFIG = {
    "type": "execution",
    "name": "pw_scf",
    "flowchartId": "unit-flowchart-id",
    "head": True,
}
DFT_METHOD_CONFIG = {"type": "pseudopotential", "subtype": "us"}
DFT_MODEL_CONFIG_WITHOUT_FUNCTIONAL = {"type": "dft", "subtype": "gga", "method": DFT_METHOD_CONFIG}


def test_creation():
    sw = Subworkflow(name=SUBWORKFLOW_NAME)
    assert sw.name == SUBWORKFLOW_NAME


@pytest.mark.parametrize("app_name", ["espresso", "vasp"])
def test_application(app_name):
    app_data = ApplicationStandata.get_by_name_first_match(app_name)
    application = Application(**app_data)
    sw = Subworkflow(name=SUBWORKFLOW_NAME, application=application)
    assert sw.application.name == app_name
    assert sw.application.version == app_data["version"]


@pytest.mark.parametrize(
    "model_type,model_subtype",
    [
        ("dft", "gga"),
        ("dft", "lda"),
    ],
)
def test_model(model_type, model_subtype):
    method = Method(type="pseudopotential", subtype="us")
    model = Model(type=model_type, subtype=model_subtype, method=method)
    sw = Subworkflow(name=SUBWORKFLOW_NAME, model=model)
    assert sw.model.type == model_type
    assert sw.model.subtype == model_subtype


@pytest.mark.parametrize("config", [DFT_MODEL_CONFIG_WITHOUT_FUNCTIONAL])
def test_model_assignment_is_coerced_to_dft_model_with_default_functional(config):
    subworkflow = Subworkflow(name=SUBWORKFLOW_NAME)
    subworkflow.model = Model(**config)
    assert isinstance(subworkflow.model, DFTModel)
    assert subworkflow.model.to_dict().get("functional") == "pbe"


def test_with_units():
    unit = Unit(**UNIT_CONFIG)
    sw = Subworkflow(name=SUBWORKFLOW_NAME, units=[unit])
    assert len(sw.units) == 1
    assert sw.units[0].name == UNIT_CONFIG["name"]


def test_id_generation():
    sw1 = Subworkflow(name=SUBWORKFLOW_NAME)
    sw2 = Subworkflow(name=SUBWORKFLOW_NAME)
    assert sw1.id != sw2.id


def test_get_as_unit():
    sw = Subworkflow(name=SUBWORKFLOW_NAME)
    unit = sw.get_as_unit()
    assert unit.type == "subworkflow"
    assert unit.id == sw.id
    assert unit.to_dict().get("_id") == sw.id
    assert unit.name == sw.name


@pytest.mark.parametrize("method", ["only_new_unit", "with_unit_instance", "with_flowchart_id"])
def test_set_unit_keeps_rendered_input_for_context_only_update(method):
    workflows = WORKFLOW_STANDATA.get_by_categories(APPLICATION_ESPRESSO, DEFAULT_WF_NAME)
    if not workflows:
        pytest.skip(f"No {DEFAULT_WF_NAME} workflow found for {APPLICATION_ESPRESSO}")

    workflow_config = workflows[0]
    wf = Workflow(**workflow_config)
    wf.add_relaxation()

    relaxation_subworkflow = wf.subworkflows[0]
    unit_to_modify = relaxation_subworkflow.get_unit_by_name(name_regex="relax")
    assert unit_to_modify is not None
    assert isinstance(unit_to_modify, ExecutionUnit)

    original_rendered = unit_to_modify.input[0].rendered

    unit_to_modify.add_context({"name": "test_key", "data": "test_value"})
    unit_to_modify.add_context({"name": "another_key", "data": 42})
    points_grid_provider = PointsGridDataProvider(
        dimensions=[2, 2, 1],
        material=Material.create(Materials.get_by_name_first_match("Silicon")),
        isEdited=True,
    )
    unit_to_modify.add_context(points_grid_provider.get_context_item_data())

    if method == "only_new_unit":
        success = relaxation_subworkflow.set_unit(unit_to_modify)
    elif method == "with_unit_instance":
        original_unit = relaxation_subworkflow.get_unit_by_name(name_regex="relax")
        success = relaxation_subworkflow.set_unit(unit_to_modify, unit=original_unit)
    elif method == "with_flowchart_id":
        flowchart_id = unit_to_modify.flowchartId
        success = relaxation_subworkflow.set_unit(unit_to_modify, unit_flowchart_id=flowchart_id)

    assert success is True

    updated_unit = relaxation_subworkflow.get_unit_by_name(name_regex="relax")
    assert updated_unit.get_context_item_data("test_key") == "test_value"
    assert updated_unit.get_context_item_data("another_key") == 42
    assert updated_unit.get_context_item_data("kgrid")["dimensions"] == [2, 2, 1]
    assert updated_unit.input[0].rendered == original_rendered
