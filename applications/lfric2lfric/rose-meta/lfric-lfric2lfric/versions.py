import re
import sys

from metomi.rose.upgrade import MacroUpgrade


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro
class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>
    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"
    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""


class vn20_t85(MacroUpgrade):
    """Upgrade macro for ticket #85 by Chris Smith."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t85"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/gungho/rose-meta/lfric-gungho
        """Add new geostrophic_forcing namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:formulation",
            r"namelist:formulation" + "\n" + " (namelist:geostrophic_forcing)",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add geostrophic_forcing setting to external_forcing namelist"""
        self.add_setting(
            config,
            ["namelist:external_forcing", "geostrophic_forcing"],
            ".false.",
        )
        """Add default data for geostrophic_forcing namelist"""
        self.add_setting(config, ["namelist:geostrophic_forcing"])
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "coordinate"], "'height'"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "heights"], "0.0"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "number_heights"], "1"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "number_times"], "1"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "profile_data_u"], "0.0"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "profile_data_v"], "0.0"
        )
        self.add_setting(
            config, ["namelist:geostrophic_forcing", "times"], "0.0"
        )
        return config, self.reports


class vn20_t358(MacroUpgrade):
    """Upgrade macro for ticket #358 by Joshua Dendy."""

    BEFORE_TAG = "vn2.0_t85"
    AFTER_TAG = "vn2.0_t358"

    def upgrade(self, config, meta_config=None):
        # Commands From: components/driver/rose-meta/lfric-driver
        """
        Add element_order_h and element_order_v to namelist finite_element,
        replacing element_order
        """
        self.add_setting(
            config, ["namelist:finite_element", "element_order_h"], "0"
        )
        self.add_setting(
            config, ["namelist:finite_element", "element_order_v"], "0"
        )
        self.remove_setting(
            config, ["namelist:finite_element", "element_order"]
        )
        return config, self.reports


class vn20_t82(MacroUpgrade):
    """Upgrade macro for ticket #82 by Chris Smith."""

    BEFORE_TAG = "vn2.0_t358"
    AFTER_TAG = "vn2.0_t82"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/gungho/rose-meta/lfric-gungho
        # Commands From: science/gungho/rose-meta/lfric-gungho
        """Add vapour_forcing namelist to configuration source list"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = source + "\n" + " (namelist:vapour_forcing)"
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        """Add vapour_forcing setting to external_forcing namelist"""
        self.add_setting(
            config, ["namelist:external_forcing", "vapour_forcing"], "'none'"
        )
        """Data for vapour_forcing namelist"""
        self.add_setting(config, ["namelist:vapour_forcing"])
        self.add_setting(
            config, ["namelist:vapour_forcing", "coordinate"], "'height'"
        )
        self.add_setting(config, ["namelist:vapour_forcing", "heights"], "0.0")
        self.add_setting(
            config, ["namelist:vapour_forcing", "number_heights"], "1"
        )
        self.add_setting(
            config, ["namelist:vapour_forcing", "number_times"], "1"
        )
        self.add_setting(
            config, ["namelist:vapour_forcing", "profile_data"], "0.0"
        )
        self.add_setting(config, ["namelist:vapour_forcing", "times"], "0.0")
        return config, self.reports


class vn20_t467(MacroUpgrade):
    """Upgrade macro for ticket #467 by Mike Whitall."""

    BEFORE_TAG = "vn2.0_t82"
    AFTER_TAG = "vn2.0_t467"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/um_physics_interface/rose-meta/um-microphysics
        # Set l_mcr_precfrac (switch for prognostic precip fraction);
        # So-far this has been hardwired in the code to be true if using
        # the comorph convection scheme, and false if not.
        # Replicate this when upgrading existing workflows from vn2.0,
        # to ensure they preserve answers
        nml = "namelist:convection"
        cv_scheme = self.get_setting_value(config, [nml, "cv_scheme"])
        if cv_scheme == "'comorph'":
            l_mcr_precfrac = ".true."
        else:
            l_mcr_precfrac = ".false."
        # end if
        # Add new microphysics settings to namelist
        nml = "namelist:microphysics"
        self.add_setting(config, [nml, "i_update_precfrac"], "'homog'")
        self.add_setting(config, [nml, "l_mcr_precfrac"], l_mcr_precfrac)
        self.add_setting(config, [nml, "l_proc_fluxes"], ".false.")
        return config, self.reports


class vn20_t339(MacroUpgrade):
    """Upgrade macro for ticket #339 by Joshua Dendy."""

    BEFORE_TAG = "vn2.0_t467"
    AFTER_TAG = "vn2.0_t339"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/gungho/rose-meta/lfric-gungho
        self.add_setting(
            config, ["namelist:idealised", "perturb_init"], ".false."
        )
        self.add_setting(
            config, ["namelist:idealised", "perturb_magnitude"], "0"
        )
        self.add_setting(config, ["namelist:idealised", "perturb_seed"], "0")
        return config, self.reports


class vn20_t541(MacroUpgrade):
    """Upgrade macro for ticket #541 by James Manners."""

    BEFORE_TAG = "vn2.0_t339"
    AFTER_TAG = "vn2.0_t541"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/socrates_interface/rose-meta/socrates-radiation
        self.add_setting(
            config, ["namelist:radiation", "cloud_entrapment"], "'zero'"
        )
        return config, self.reports


class vn20_t481(MacroUpgrade):
    """Upgrade macro for ticket #481 by Mike Whitall."""

    BEFORE_TAG = "vn2.0_t541"
    AFTER_TAG = "vn2.0_t481"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/um_physics_interface/rose-meta/um-cloud
        # All altered settings are in the cloud namelist
        nml = "namelist:cloud"
        # Logical ez_subcrit_only replaced by a 3-way integer switch:
        ez_subcrit = self.get_setting_value(config, [nml, "ez_subcrit"])
        self.remove_setting(config, [nml, "ez_subcrit"])
        if ez_subcrit == ".true.":
            i_bm_ez = "'subcrit'"
        else:
            i_bm_ez = "'orig'"
        # end if
        self.add_setting(config, [nml, "i_bm_ez_opt"], i_bm_ez)
        # New settings added
        self.add_setting(config, [nml, "ent_coef_bm"], "0.2")
        self.add_setting(config, [nml, "l_bm_sigma_s_grad"], ".false.")
        self.add_setting(config, [nml, "l_bm_tweaks"], ".false.")
        self.add_setting(config, [nml, "max_sigmas"], "3.0")
        self.add_setting(config, [nml, "min_sigx_ft"], "0.0")
        self.add_setting(config, [nml, "turb_var_fac_bm"], "1.0")
        return config, self.reports


class vn20_t450(MacroUpgrade):
    """Upgrade macro for ticket #450 by Juan M Castillo."""

    BEFORE_TAG = "vn2.0_t481"
    AFTER_TAG = "vn2.0_t450"

    def upgrade(self, config, meta_config=None):
        # Commands From: applications/lfric2lfric/rose-meta/lfric-lfric2lfric
        # Commands From: applications/lfric2lfric/rose-meta/lfric-lfric2lfric
        # Commands From: applications/lfric2lfric/rose-meta/lfric-lfric2lfric
        self.add_setting(
            config,
            ["namelist:lfric2lfric", "destination_mesh_name"],
            "'dynamics'",
        )
        self.add_setting(
            config,
            ["namelist:lfric2lfric", "destination_meshfile_prefix"],
            "'input_mesh_faces'",
        )
        self.add_setting(
            config,
            ["namelist:lfric2lfric", "source_mesh_name"],
            "'multigrid_l1'",
        )
        self.add_setting(
            config,
            ["namelist:lfric2lfric", "source_meshfile_prefix"],
            "'input_mesh_faces'",
        )
        self.remove_setting(config, ["namelist:lfric2lfric", "chain_mesh_tags"])

        return config, self.reports
