import { Template } from "@exabyte-io/ade.js";

import { wodeProviders } from "./context/providers";

// We patch the static providerRegistry here so that
// Template has all context providers available
// to it when creating workflows. It is then re-exported
// from WoDe for use downstream.

Template.setContextProvidersConfig(wodeProviders);

export { Template };
